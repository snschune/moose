//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RayTracingStudy.h"

// Local Includes
#include "AuxRayKernel.h"
#include "RayBoundaryConditionBase.h"
#include "RayKernel.h"
#include "TraceRay.h"
#include "TraceRayTools.h"

// MOOSE Includes
#include "AuxiliarySystem.h"
#include "Assembly.h"
#include "NonlinearSystemBase.h"
#include "TimedPrint.h"
#include "LIFOBuffer.h"
#include "CircularBuffer.h"

// libMesh Includes
#include "libmesh/enum_to_string.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/parallel_sync.h"

InputParameters
RayTracingStudy::validParams()
{
  auto params = GeneralUserObject::validParams();

  // Parameters for the execution of the rays
  params += ParallelRayStudy::validParams();

  params.addRangeCheckedParam<Real>("ray_distance",
                                    std::numeric_limits<Real>::max(),
                                    "ray_distance > 0",
                                    "The maximum distance a Ray can travel");

  params.addParam<bool>(
      "tolerate_failure", false, "Whether or not to tolerate a ray tracing failure");

  MooseEnum work_buffers("lifo circular", "circular");
  params.addParam<MooseEnum>("work_buffer_type", work_buffers, "The work buffer type to use");

  params.addParam<bool>(
      "ray_kernel_coverage_check", true, "Whether or not to perform coverage checks on RayKernels");
  params.addParam<bool>("planar_face_check",
                        true,
                        "Whether or not to check and warn if any element's faces are not planar");

  params.addParam<bool>(
      "always_cache_traces",
      false,
      "Whether or not to cache the Ray traces on every execution, primarily for use in output. "
      "Warning: this can get expensive very quick with a large number of rays!");
  params.addParam<bool>("data_on_cache_traces",
                        false,
                        "Whether or not to also cache the Ray's data when caching its traces");
  params.addParam<bool>("aux_data_on_cache_traces",
                        false,
                        "Whether or not to also cache the Ray's aux data when caching its traces");
  params.addParam<bool>(
      "segments_on_cache_traces",
      true,
      "Whether or not to cache individual segments when trace caching is enabled. If false, we "
      "will instead cache a segment for each part of the trace where the direction is the same. "
      "This minimizes the number of segments requied to represent the Ray's path, but removes the "
      "ability to show Ray field data on each segment through an element.");

  params.addParam<bool>("use_internal_sidesets",
                        false,
                        "Whether or not to use internal sidesets for RayBCs in ray tracing");

  params.addParam<bool>("warn_subdomain_hmax",
                        true,
                        "Whether or not to warn if the approximated hmax (constant on subdomain) "
                        "varies significantly for an element");

  params.addParam<bool>(
      "verify_rays",
      true,
      "Whether or not to verify if Rays have valid information in optmized modes before being "
      "traced. All Rays are verified outside of optimized modes regardless of this parameter.");

  ExecFlagEnum & exec_enum = params.set<ExecFlagEnum>("execute_on", true);
  exec_enum.addAvailableFlags(EXEC_PRE_KERNELS);

  // Whether or not each Ray must be registered using the registerRay() API
  params.addPrivateParam<bool>("_use_ray_registration", true);
  // Whether or not to bank Rays on completion
  params.addPrivateParam<bool>("_bank_rays_on_completion", true);
  /// Whether or not subdomain setup is dependent on the Ray
  params.addPrivateParam<bool>("_ray_dependent_subdomain_setup", true);

  // Add a point neighbor relationship manager
  params.addRelationshipManager("ElementPointNeighborLayers",
                                Moose::RelationshipManagerType::GEOMETRIC |
                                    Moose::RelationshipManagerType::ALGEBRAIC,
                                [](const InputParameters &, InputParameters & rm_params) {
                                  rm_params.set<unsigned short>("layers") = 1;
                                });

  return params;
}

RayTracingStudy::RayTracingStudy(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    SidePtrHelper(),
    _mesh(_fe_problem.mesh()),
    _comm(_mesh.comm()),
    _pid(_comm.rank()),
    _error_prefix(type() + " '" + name() + "'"),

    _ray_kernel_coverage_check(getParam<bool>("ray_kernel_coverage_check")),
    _planar_face_check(getParam<bool>("planar_face_check")),
    _use_ray_registration(getParam<bool>("_use_ray_registration")),
    _use_internal_sidesets(getParam<bool>("use_internal_sidesets")),
    _tolerate_failure(getParam<bool>("tolerate_failure")),
    _bank_rays_on_completion(getParam<bool>("_bank_rays_on_completion")),
    _ray_dependent_subdomain_setup(getParam<bool>("_ray_dependent_subdomain_setup")),

    _always_cache_traces(getParam<bool>("always_cache_traces")),
    _data_on_cache_traces(getParam<bool>("data_on_cache_traces")),
    _aux_data_on_cache_traces(getParam<bool>("aux_data_on_cache_traces")),
    _segments_on_cache_traces(getParam<bool>("segments_on_cache_traces")),
    _ray_max_distance(getParam<Real>("ray_distance")),
    _verify_rays(getParam<bool>("verify_rays")),

    _execute_study_timer(registerTimedSection("executeStudy", 1)),
    _generate_timer(registerTimedSection("generate", 1)),
    _propagate_timer(registerTimedSection("propagate", 1)),

    _cached_traces(nullptr),
    _threaded_cached_traces(libMesh::n_threads()),

    _computing_jacobian(false),
    _computing_residual(false),
    _num_cached(libMesh::n_threads(), 0),

    _has_internal_sidesets(false),
    _internal_sideset_elem_integer_name(name() + "_internal_sideset"),
    _has_same_level_active_elems(sameLevelActiveElems()),

    _b_box(MeshTools::create_nodal_bounding_box(_mesh.getMesh())),
    _domain_max_length(1.01 * (_b_box.max() - _b_box.min()).norm()),
    _total_volume(computeTotalVolume()),

    _threaded_cache_ray_kernel(libMesh::n_threads()),
    _threaded_cache_ray_bc(libMesh::n_threads()),
    _threaded_ray_object_registration(libMesh::n_threads()),
    _threaded_current_ray_kernels(libMesh::n_threads()),
    _threaded_trace_ray(libMesh::n_threads()),
    _threaded_fe_face(libMesh::n_threads()),
    _threaded_q_face(libMesh::n_threads()),
    _threaded_cached_normals(libMesh::n_threads()),
    _threaded_next_ray_id(libMesh::n_threads()),

    _parallel_ray_study(libmesh_make_unique<ParallelRayStudy>(*this, _threaded_trace_ray)),

    _local_trace_ray_results(TraceRay::FAILED_TRACES + 1, 0),

    _called_initial_setup(false)
{
  // Initialize a tracing object for each thread
  for (THREAD_ID tid = 0; tid < libMesh::n_threads(); ++tid)
  {
    // Initialize a tracing object for each thread
    _threaded_trace_ray[tid] = std::make_shared<TraceRay>(*this, tid);

    // Setup the face FEs for normal computation on the fly
    _threaded_fe_face[tid] = FEBase::build(_mesh.dimension(), FEType(CONSTANT, MONOMIAL));
    _threaded_q_face[tid] = QBase::build(QGAUSS, _mesh.dimension() - 1, CONSTANT);
    _threaded_fe_face[tid]->attach_quadrature_rule(_threaded_q_face[tid].get());
    _threaded_fe_face[tid]->get_normals();
  }

  if (_execute_enum.size() > 1 && _execute_enum.contains(EXEC_PRE_KERNELS))
    paramError(
        "execute_on",
        "PRE_KERNELS cannot be mixed with any other execution flag.\nThat is, you cannot currently "
        "mix RayKernels that contribute to the Jacobian/residual with those that do not.");

  resetUniqueRayIDs();
}

void
RayTracingStudy::initialSetup()
{
  // Keep track of initialSetup call to avoid registration of various things
  _called_initial_setup = true;

  // Check for RayKernel coverage
  coverageChecks();

  // Make sure the dependencies exist, if any
  dependencyChecks();

  // Check for traceable element types
  traceableMeshChecks();

  // Setup Ray -> RayTracingObject maps if we have registered Rays
  associateRegisteredRays();

  // Setup for internal sidesets
  internalSidesetSetup();

  // Setup approximate hmax for each subdomain
  setupSubdomainHmax();

  // Call initial setup on all of the objects
  for (auto & rto : getRayTracingObjects())
    rto->initialSetup();

  // Check for proper exec flags with RayKernels
  std::vector<RayKernelBase *> ray_kernels;
  getRayKernels(ray_kernels, 0);
  for (const auto & rkb : ray_kernels)
    if (dynamic_cast<RayKernel *>(rkb) && !_execute_enum.contains(EXEC_PRE_KERNELS))
      mooseError(_error_prefix,
                 ": Has RayKernel objects that contribute to residuals and Jacobians.",
                 "\nIn this case, the study should use the execute_on = PRE_KERNELS");

  // Build 1D quadrature rule for along a segment
  _segment_qrule =
      QBase::build(QGAUSS, 1, _fe_problem.getNonlinearSystemBase().getMinQuadratureOrder());
}

void
RayTracingStudy::residualSetup()
{
  _computing_residual = true;
  mooseAssert(_num_cached == std::vector<unsigned int>(libMesh::n_threads(), 0),
              "Should be all zero");

  for (auto & rto : getRayTracingObjects())
    rto->residualSetup();
}

void
RayTracingStudy::jacobianSetup()
{
  _computing_jacobian = true;
  mooseAssert(_num_cached == std::vector<unsigned int>(libMesh::n_threads(), 0),
              "Should be all zero");

  for (auto & rto : getRayTracingObjects())
    rto->jacobianSetup();
}

void
RayTracingStudy::timestepSetup()
{
  for (auto & rto : getRayTracingObjects())
    rto->timestepSetup();
}

void
RayTracingStudy::meshChanged()
{
  // Internal sidesets may have moved
  internalSidesetSetup();

  setupSubdomainHmax();

  _has_same_level_active_elems = sameLevelActiveElems();

  for (const auto & trace_ray : _threaded_trace_ray)
    trace_ray->meshChanged();
}

void
RayTracingStudy::execute()
{
  executeStudy();
}

void
RayTracingStudy::coverageChecks()
{
  // Check for coverage of RayKernels on domain
  if (_ray_kernel_coverage_check)
  {
    std::vector<RayKernelBase *> ray_kernels;
    getRayKernels(ray_kernels, 0);

    std::set<SubdomainID> ray_kernel_blocks;
    for (const auto & rk : ray_kernels)
      ray_kernel_blocks.insert(rk->blockIDs().begin(), rk->blockIDs().end());

    std::set<SubdomainID> missing;
    std::set_difference(_mesh.meshSubdomains().begin(),
                        _mesh.meshSubdomains().end(),
                        ray_kernel_blocks.begin(),
                        ray_kernel_blocks.end(),
                        std::inserter(missing, missing.begin()));

    if (!missing.empty() && !ray_kernel_blocks.count(Moose::ANY_BLOCK_ID))
    {
      std::ostringstream error;
      error << _error_prefix << ": Subdomains { ";
      std::copy(missing.begin(), missing.end(), std::ostream_iterator<SubdomainID>(error, " "));
      error << "} do not have RayKernels defined!";

      mooseError(error.str());
    }
  }
}

void
RayTracingStudy::dependencyChecks()
{
  std::vector<RayTracingObject *> ray_tracing_objects;

  getRayKernels(ray_tracing_objects, 0);
  verifyDependenciesExist(ray_tracing_objects);

  getRayBCs(ray_tracing_objects, 0);
  verifyDependenciesExist(ray_tracing_objects);
}

void
RayTracingStudy::verifyDependenciesExist(const std::vector<RayTracingObject *> & rtos)
{
  for (const auto & rto : rtos)
    for (const auto & dep_name : rto->getRequestedItems())
    {
      bool found = false;
      for (const auto & rto_search : rtos)
        if (rto_search->name() == dep_name)
        {
          found = true;
          break;
        }

      if (!found)
      {
        const auto object_type = rto->parameters().get<std::string>("_moose_warehouse_system_name");

        mooseError(rto->type(),
                   " '",
                   rto->name(),
                   "' depends on ",
                   object_type,
                   " '",
                   dep_name,
                   "' but such an object does not exist");
      }
    }
}

void
RayTracingStudy::traceableMeshChecks()
{
  bool planar_face_check_warned = false;

  for (const auto & elem : *_mesh.getActiveLocalElementRange())
  {
    if (_fe_problem.adaptivity().isOn())
    {
      if (!TraceRayTools::isAdaptivityTraceableElem(elem))
        mooseError(_error_prefix,
                   ": Element type ",
                   Utility::enum_to_string(elem->type()),
                   " is not supported in ray tracing with adaptivity");
    }
    else
    {
      if (!TraceRayTools::isTraceableElem(elem))
        mooseError(_error_prefix,
                   ": Element type ",
                   Utility::enum_to_string(elem->type()),
                   " is not supported in ray tracing");
    }

    if (_planar_face_check && !planar_face_check_warned)
      for (const auto s : elem->side_index_range())
        if (!sidePtrHelper(elem, s)->has_affine_map())
        {
          planar_face_check_warned = true;
          mooseWarning(_error_prefix,
                       ":\n\nThe mesh contains element(s) that do not have planar faces.",
                       "\nIf the faces are sufficiently non-planar, ray tracing may fail.",
                       "\n\nTo disable this warning, set planar_face_check = false.");
          break;
        }
  }
}

void
RayTracingStudy::internalSidesetSetup()
{
  // Even if we have _use_internal_sidesets == false, we will make sure the user didn't add RayBCs
  // on internal boundaries
  _has_internal_sidesets = false;
  _internal_sidesets.clear();

  // First, we are going to store all elements with internal sidesets (if any) that have active
  // RayBCs on them as elem -> vector of (side, vector of boundary ids)
  std::unordered_map<Elem *, std::vector<std::pair<unsigned short, std::vector<BoundaryID>>>> map;
  for (const auto & bnd_elem : *_mesh.getBoundaryElementRange())
  {
    Elem * elem = bnd_elem->_elem;
    const unsigned int side = bnd_elem->_side;
    const auto bnd_id = bnd_elem->_bnd_id;

    // Not internal
    const Elem * const neighbor = elem->neighbor_ptr(side);
    if (!neighbor)
      continue;

    // No RayBCs on this sideset
    std::vector<RayBoundaryConditionBase *> result;
    getRayBCs(result, bnd_id, 0);
    if (result.empty())
      continue;

    if (neighbor->subdomain_id() == elem->subdomain_id())
      mooseError(_error_prefix,
                 ":\n\nRayBCs exist on internal sidesets that are not bounded by a different",
                 "\nsubdomain on each side.",
                 "\n\nIn order to use RayBCs on internal sidesets, said sidesets must have",
                 "\na different subdomain on each side.");

    auto & map_entry = map[elem];
    bool found_side = false;
    for (auto & pair : map_entry)
      if (pair.first == side)
      {
        pair.second.push_back(bnd_id);
        found_side = true;
      }
    if (!found_side)
      map_entry.emplace_back(side, std::vector<BoundaryID>(1, bnd_id));
  }

  // Now, we need to determine if anyone has internal sidesets that also have RayBCs. If so, we need
  // to add the extra element integer that will let us quickly index into the final storage for
  // internal sidesets. If said extra element integer is DofObject::invalid_id (the default), it
  // means that said element does not have an entry in the internal sidesets vector and therefore
  // has no internal sidesets.
  unsigned int has_internal_sidesets = !map.empty();
  _communicator.sum(has_internal_sidesets);
  if (has_internal_sidesets)
  {
    if (!_use_internal_sidesets)
      mooseError(_error_prefix,
                 ":\n\n RayBCs are defined on internal sidesets, but the study is not set to use ",
                 "internal sidesets during tracing.",
                 "\n\nSet the parameter use_internal_sidesets = true to enable this capability.");

    _has_internal_sidesets = true;

    // We added the integer previously: invalidate it everywhere
    if (_mesh.getMesh().has_elem_integer(_internal_sideset_elem_integer_name))
      for (auto elem : _mesh.getMesh().element_ptr_range())
        elem->set_extra_integer(_internal_sideset_elem_integer, DofObject::invalid_id);
    // Didn't add the elem integer previously: need to add it (it initializes as invalid)
    else
      _internal_sideset_elem_integer =
          _mesh.getMesh().add_elem_integer(_internal_sideset_elem_integer_name);

    // Move the data from the map into a continguous vector that is indexed using the extra integer
    _internal_sidesets.reserve(map.size());
    for (const auto & pair : map)
    {
      Elem * elem = pair.first;
      auto & entry = pair.second;
      _internal_sidesets.push_back(std::move(entry));
      elem->set_extra_integer(_internal_sideset_elem_integer, _internal_sidesets.size() - 1);
    }
  }
}

void
RayTracingStudy::setupSubdomainHmax()
{
  // Setup map with subdomain keys
  _subdomain_hmax.clear();
  for (const auto subdomain_id : _mesh.meshSubdomains())
    _subdomain_hmax[subdomain_id] = std::numeric_limits<Real>::min();

  // Set local max for each subdomain
  for (const auto & elem : *_mesh.getActiveLocalElementRange())
  {
    auto & entry = _subdomain_hmax.at(elem->subdomain_id());
    entry = std::max(entry, elem->hmax());
  }

  // Accumulate global max for each subdomain
  _communicator.max(_subdomain_hmax);

  if (getParam<bool>("warn_subdomain_hmax"))
  {
    const auto warn_prefix = type() + " '" + name() + "': ";
    const auto warn_suffix =
        "\n\nRay tracing uses an approximate element size for each subdomain to scale the\n"
        "tolerances used in computing ray intersections. This warning suggests that the\n"
        "approximate element size is not a good approximation. This is likely due to poor\n"
        "element aspect ratios.\n\n"
        "To disable this warning, set warn_subdomain_hmax = false.\n";

    for (const auto & elem : *_mesh.getActiveLocalElementRange())
    {
      const auto hmin = elem->hmin();
      const auto hmax = elem->hmax();
      const auto max_hmax = subdomainHmax(elem->subdomain_id());

      const auto hmax_rel = hmax / max_hmax;
      if (hmax_rel < 1.e-2 || hmax_rel > 1.e2)
        mooseWarning(
            warn_prefix, "Element hmax varies significantly from subdomain hmax.\n", warn_suffix);

      const auto h_rel = max_hmax / hmin;
      if (h_rel > 1.e2)
        mooseWarning(
            warn_prefix, "Element hmin varies significantly from subdomain hmax.\n", warn_suffix);
    }
  }
}

void
RayTracingStudy::associateRegisteredRays()
{
  // First, clear the objects associated with each Ray on each thread
  for (THREAD_ID tid = 0; tid < libMesh::n_threads(); ++tid)
    for (auto & set : _threaded_ray_object_registration[tid])
      set.clear();

  for (auto & rto : getRayTracingObjects())
  {
    const auto object_type = rto->parameters().get<std::string>("_moose_warehouse_system_name");
    const auto & params = rto->parameters();

    // The Ray names associated with this RayTracingObject
    const auto & ray_names = params.get<std::vector<std::string>>("rays");
    // The registration for RayTracingObjects for the thread rto is on
    const auto tid = params.get<THREAD_ID>("_tid");
    auto & registration = _threaded_ray_object_registration[tid];

    // No ray names but we need them
    if (ray_names.empty() && _use_ray_registration)
      mooseError(_error_prefix,
                 ": ",
                 object_type,
                 " '",
                 rto->name(),
                 "' must supply the Rays that it is associated with via the rays parameter.\n",
                 "This requirement is set by the '_use_ray_registration' private parameter.");
    // Ray names but we don't need them
    if (!ray_names.empty() && !_use_ray_registration)
      mooseError(_error_prefix,
                 ": ",
                 object_type,
                 " '",
                 rto->name(),
                 "' has supplied the Rays that it is associated with, but the study does not "
                 "require Ray registration");

    // Register each Ray for this object in the registration
    for (const auto & ray_name : ray_names)
    {
      const auto id = registeredRayID(ray_name, /* graceful = */ true);
      if (id == DofObject::invalid_id)
        mooseError(_error_prefix,
                   ": The Ray '",
                   ray_name,
                   "' requested by ",
                   object_type,
                   " '",
                   rto->name(),
                   "' is not a registered Ray");
      registration[id].insert(rto);
    }
  }
}

void
RayTracingStudy::zeroAuxVariables()
{
  // TODO: If we have Rays with different execute_on flags, make sure to zero only the relavent
  // aux vars
  std::set<std::string> vars_to_be_zeroed;

  std::vector<RayKernelBase *> ray_kernels;
  getRayKernels(ray_kernels, 0);
  for (auto & rk : ray_kernels)
  {
    AuxRayKernel * aux_rk = dynamic_cast<AuxRayKernel *>(rk);
    if (aux_rk)
      vars_to_be_zeroed.insert(aux_rk->variable().name());
  }

  std::vector<std::string> vars_to_be_zeroed_vec(vars_to_be_zeroed.begin(),
                                                 vars_to_be_zeroed.end());
  _fe_problem.getAuxiliarySystem().zeroVariables(vars_to_be_zeroed_vec);
}

void
RayTracingStudy::subdomainSetup(SubdomainID subdomain, THREAD_ID tid, RayID ray_id)
{
  // Call subdomain setup on FE
  _fe_problem.subdomainSetup(subdomain, tid);

  std::set<MooseVariableFEBase *> needed_moose_vars;
  std::set<unsigned int> needed_mat_props;

  // Get RayKernels and their dependencies and call subdomain setup
  getRayKernels(_threaded_current_ray_kernels[tid], subdomain, tid, ray_id);
  for (auto & rkb : _threaded_current_ray_kernels[tid])
  {
    rkb->subdomainSetup();

    const auto & mv_deps = rkb->getMooseVariableDependencies();
    needed_moose_vars.insert(mv_deps.begin(), mv_deps.end());

    const auto & mp_deps = rkb->getMatPropDependencies();
    needed_mat_props.insert(mp_deps.begin(), mp_deps.end());
  }

  // Prepare aux vars
  for (auto & var : needed_moose_vars)
    if (var->kind() == Moose::VarKindType::VAR_AUXILIARY)
      var->prepareAux();

  _fe_problem.setActiveElementalMooseVariables(needed_moose_vars, tid);
  _fe_problem.setActiveMaterialProperties(needed_mat_props, tid);
  _fe_problem.prepareMaterials(subdomain, tid);
}

void
RayTracingStudy::reinitSegment(
    const Elem * elem, const Point & start, const Point & end, const Real length, THREAD_ID tid)
{
  mooseAssert(MooseUtils::absoluteFuzzyEqual((start - end).norm(), length), "Invalid length");

  _fe_problem.setCurrentSubdomainID(elem, tid);

  // If we have any variables or material properties that are active, we definitely need to reinit
  bool reinit = _fe_problem.hasActiveElementalMooseVariables(tid) ||
                _fe_problem.hasActiveMaterialProperties(tid);
  // If not, make sure that the RayKernels have not requested a reinit (this could happen when a
  // RayKernel doesn't have variables or materials but still does an integration and needs qps)
  if (!reinit)
    for (const RayKernelBase * rk : currentRayKernels(tid))
      if (rk->needSegmentReinit())
      {
        reinit = true;
        break;
      }

  if (reinit)
  {
    _fe_problem.prepare(elem, tid);

    std::vector<Point> points;
    std::vector<Real> weights;
    buildSegmentQuadrature(start, end, length, points, weights);
    _fe_problem.reinitElemPhys(elem, points, tid);
    _fe_problem.assembly(tid).modifyArbitraryWeights(weights);

    _fe_problem.reinitMaterials(elem->subdomain_id(), tid);
  }
}

void
RayTracingStudy::buildSegmentQuadrature(const Point & start,
                                        const Point & end,
                                        const Real length,
                                        std::vector<Point> & points,
                                        std::vector<Real> & weights) const
{
  points.resize(_segment_qrule->n_points());
  weights.resize(_segment_qrule->n_points());

  const Point diff = end - start;
  const Point sum = end + start;
  mooseAssert(MooseUtils::absoluteFuzzyEqual(length, diff.norm()), "Invalid length");

  // The standard quadrature rule should be on x = [-1, 1]
  // To scale the points, you...
  //  - Scale to size of the segment in 3D
  //    initial_scaled_qp = x_qp * 0.5 * (end - start) = 0.5 * x_qp * diff
  //  - Shift quadrature midpoint to segment midpoint
  //    final_qp = initial_scaled_qp + 0.5 * (end - start) = initial_scaled_qp + 0.5 * sum
  //             = 0.5 * (x_qp * diff + sum)
  for (unsigned int qp = 0; qp < _segment_qrule->n_points(); ++qp)
  {
    points[qp] = 0.5 * (_segment_qrule->qp(qp)(0) * diff + sum);
    weights[qp] = 0.5 * _segment_qrule->w(qp) * length;
  }
}

void
RayTracingStudy::postOnSegment(THREAD_ID tid)
{
  mooseAssert(_num_cached[tid] != 0 ? _computing_jacobian || _computing_residual : true,
              "Values should only be cached when computing Jacobian/residual");

  // Fill into cached Jacobian/residuals if necessary
  if (_computing_jacobian)
  {
    _fe_problem.cacheJacobian(tid);

    if (++_num_cached[tid] == 20)
    {
      Threads::spin_mutex::scoped_lock lock(_spin_mutex);
      _fe_problem.addCachedJacobian(tid);
      _num_cached[tid] = 0;
    }
  }
  else if (_computing_residual)
  {
    _fe_problem.cacheResidual(tid);

    if (++_num_cached[tid] == 20)
    {
      Threads::spin_mutex::scoped_lock lock(_spin_mutex);
      _fe_problem.addCachedResidual(tid);
      _num_cached[tid] = 0;
    }
  }
}

void
RayTracingStudy::preExecuteStudy()
{
  // Zero the AuxVariables that our AuxRayKernels contribute to before they accumulate
  zeroAuxVariables();

  for (auto & rto : getRayTracingObjects())
    rto->preExecuteStudy();
}

void
RayTracingStudy::postExecuteStudy()
{
  // Add any stagglers that contribute to the Jacobian or residual
  for (THREAD_ID tid = 0; tid < libMesh::n_threads(); ++tid)
    if (_num_cached[tid] != 0)
    {
      if (_computing_jacobian)
        _fe_problem.addCachedJacobian(tid);
      else if (_computing_residual)
        _fe_problem.addCachedResidual(tid);
      else
        mooseError("Should not have cached values without Jacobian/residual computation");

      _num_cached[tid] = 0;
    }

  _computing_jacobian = false;
  _computing_residual = false;

  // AuxRayKernels may have modified AuxVariables
  _fe_problem.getAuxiliarySystem().solution().close();
  _fe_problem.getAuxiliarySystem().update();

  // Clear FE
  for (THREAD_ID tid = 0; tid < libMesh::n_threads(); ++tid)
  {
    _fe_problem.clearActiveElementalMooseVariables(tid);
    _fe_problem.clearActiveMaterialProperties(tid);
  }

  for (auto & rto : getRayTracingObjects())
    rto->postExecuteStudy();
}

void
RayTracingStudy::executeStudy()
{
  TIME_SECTION(_execute_study_timer);

  if (!_called_initial_setup)
    mooseError(_error_prefix,
               " did not call RayTracingStudy::initialSetup().\n\n",
               "You likely override initialSetup() and forgot to call the parent's method.");

  // Reset ray start/complete timers
  _total_intersections = 0;
  _max_intersections = 0;
  _max_trajectory_changes = 0;

  // Reset physical tracing stats
  for (auto & val : _local_trace_ray_results)
    val = 0;

  // Reset crossing and intersection
  _ending_processor_crossings = 0;
  _total_processor_crossings = 0;
  _ending_max_processor_crossings = 0;
  _ending_intersections = 0;
  _ending_max_intersections = 0;
  _ending_max_trajectory_changes = 0;
  _ending_distance = 0;
  _total_distance = 0;

  preExecuteStudy();

  _ray_bank.clear();
#ifndef NDEBUG
  _debug_ray_bank.clear();
#endif
  for (THREAD_ID tid = 0; tid < libMesh::n_threads(); ++tid)
  {
    _threaded_trace_ray[tid]->preExecute();
    _threaded_cached_normals[tid].clear();
  }

  // This releases our use of the old cached trace map. Other objects may still use it if they
  // properly obtained its shared pointer. Otherwise, it's deleted if no one else is using it.
  _cached_traces = nullptr;

  // Create a cached trace map for each thread
  for (THREAD_ID tid = 0; tid < libMesh::n_threads(); ++tid)
  {
    mooseAssert(!_threaded_cached_traces[tid], "Cached trace map was not deleted properly");
    _threaded_cached_traces[tid] = libmesh_make_unique<std::vector<TraceData>>();
  }

  _communicator.barrier();
  _execution_start_time = std::chrono::steady_clock::now();

  _parallel_ray_study->preExecute();

  {
    {
      auto generation_start_time = std::chrono::steady_clock::now();

      TIME_SECTION(_generate_timer);

      generateRays();

      _generation_time = std::chrono::steady_clock::now() - generation_start_time;
    }

#ifndef NDEBUG
    // At this point, nobody is working so this is good time to make sure
    // Ray IDs are unique across all processors in the working buffer
    verifyUniqueRayIDs(_parallel_ray_study->workingObjectBuffer().begin(),
                       _parallel_ray_study->workingObjectBuffer().end(),
                       /* global = */ true,
                       /* error_suffix = */ "after generateRays()");

    verifyUniqueRays(_parallel_ray_study->workingObjectBuffer().begin(),
                     _parallel_ray_study->workingObjectBuffer().end(),
                     /* error_suffix = */ "after generateRays()");
#endif

    {
      CONSOLE_TIMED_PRINT("Propagating rays");

      const auto propagation_start_time = std::chrono::steady_clock::now();

      TIME_SECTION(_propagate_timer);

      _parallel_ray_study->execute();

      _propagation_time = std::chrono::steady_clock::now() - propagation_start_time;
    }
  }

  _execution_time = std::chrono::steady_clock::now() - _execution_start_time;

#ifndef NDEBUG
  // Outside of debug, _debug_ray_bank holds all of the Rays that have ended on this processor
  // We can use this as a global point to check for unique IDs for every Ray that has traced
  verifyUniqueRayIDs(_debug_ray_bank.begin(),
                     _debug_ray_bank.end(),
                     /* global = */ true,
                     /* error_suffix = */ "after tracing completed");

  verifyUniqueRays(_parallel_ray_study->workingObjectBuffer().begin(),
                   _parallel_ray_study->workingObjectBuffer().end(),
                   /* error_suffix = */ "after tracing completed");
#endif

  // Update counters from the threaded trace objects
  for (const auto & tr : _threaded_trace_ray)
    for (std::size_t i = 0; i < _local_trace_ray_results.size(); ++i)
      _local_trace_ray_results[i] += tr->results()[i];

  // Update local ending counters
  _total_processor_crossings = _ending_processor_crossings;
  _max_processor_crossings = _ending_max_processor_crossings;
  _total_intersections = _ending_intersections;
  _max_intersections = _ending_max_intersections;
  _max_trajectory_changes = _ending_max_trajectory_changes;
  _total_distance = _ending_distance;
  // ...and communicate the global values
  _communicator.sum(_total_processor_crossings);
  _communicator.max(_max_processor_crossings);
  _communicator.sum(_total_intersections);
  _communicator.max(_max_intersections);
  _communicator.max(_max_trajectory_changes);
  _communicator.sum(_total_distance);

  // Throw a warning with the number of failed (tolerated) traces
  if (_tolerate_failure)
  {
    auto failures = _local_trace_ray_results[TraceRay::FAILED_TRACES];
    _communicator.sum(failures);
    if (failures)
      mooseWarning(
          type(), " '", name(), "': ", failures, " ray-tracing failures were tolerated.\n");
  }

  // Clear the current RayKernels
  for (THREAD_ID tid = 0; tid < libMesh::n_threads(); ++tid)
    _threaded_current_ray_kernels[tid].clear();

  // Whether we have one thread or multiple, the single cache trace structure will always start
  // with the map from thread 0. We will just move ownership from the threaded vector to the main
  // map. Then, if we have more threads, we will just add to the newly moved map. This avoids
  // copying some or all of the data.
  _cached_traces = std::move(_threaded_cached_traces[0]);
  _threaded_cached_traces[0] = nullptr;

  // Copy cached trace data from thread IDs > 1
  if (libMesh::n_threads() > 1)
  {
    // Reserve space for the extra entries
    auto num_entries = _cached_traces->size();
    for (THREAD_ID tid = 1; tid < libMesh::n_threads(); ++tid)
      num_entries += _threaded_cached_traces[tid]->size();
    _cached_traces->reserve(num_entries);

    // Move the contents from each threaded map into the main map
    for (THREAD_ID tid = 1; tid < libMesh::n_threads(); ++tid)
    {
      for (const auto & entry : *_threaded_cached_traces[tid])
        _cached_traces->emplace_back(std::move(entry));

      // This will delete the threaded map by resetting the unique pointer
      _threaded_cached_traces[tid] = nullptr;
    }
  }

  postExecuteStudy();
}

void
RayTracingStudy::onCompleteRay(const std::shared_ptr<Ray> & ray)
{
  mooseAssert(currentlyPropagating(), "Should only be called during Ray propagation");

  _ending_processor_crossings += ray->processorCrossings();
  _ending_max_processor_crossings =
      std::max(_ending_max_processor_crossings, ray->processorCrossings());
  _ending_intersections += ray->intersections();
  _ending_max_intersections = std::max(_ending_max_intersections, ray->intersections());
  _ending_max_trajectory_changes =
      std::max(_ending_max_trajectory_changes, ray->trajectoryChanges());
  _ending_distance += ray->distance();

  if (_bank_rays_on_completion)
    _ray_bank.push_back(ray);

#ifndef NDEBUG
  // Outside of optimized modes, we will save all Rays that end on this processor
  // in their own bank for use in debugging
  _debug_ray_bank.emplace_back(ray);
#endif
}

RayDataIndex
RayTracingStudy::registerRayDataInternal(const std::string & name, const bool aux)
{
  Threads::spin_mutex::scoped_lock lock(_spin_mutex);

  if (_called_initial_setup)
    mooseError(
        _error_prefix, ": Cannot register Ray ", (aux ? "aux" : ""), " data after initialSetup()");

  auto & map = aux ? _ray_aux_data_map : _ray_data_map;
  const auto find = map.find(name);
  if (find != map.end())
    return find->second;

  auto & other_map = aux ? _ray_data_map : _ray_aux_data_map;
  if (other_map.find(name) != other_map.end())
    mooseError(_error_prefix,
               ": Cannot register Ray aux data with name ",
               name,
               " because Ray ",
               (aux ? "(non-aux)" : "aux"),
               " data already exists with said name.");

  // Add into the name -> index map
  map.emplace(name, map.size());

  // Add into the index -> names vector
  auto & vector = aux ? _ray_aux_data_names : _ray_data_names;
  vector.push_back(name);

  return map.size() - 1;
}

std::vector<RayDataIndex>
RayTracingStudy::registerRayDataInternal(const std::vector<std::string> & names, const bool aux)
{
  std::vector<RayDataIndex> indices(names.size());
  for (std::size_t i = 0; i < names.size(); ++i)
    indices[i] = registerRayDataInternal(names[i], aux);
  return indices;
}

RayDataIndex
RayTracingStudy::getRayDataIndexInternal(const std::string & name,
                                         const bool aux,
                                         const bool graceful) const
{
  Threads::spin_mutex::scoped_lock lock(_spin_mutex);

  const auto & map = aux ? _ray_aux_data_map : _ray_data_map;
  const auto find = map.find(name);
  if (find != map.end())
    return find->second;

  if (graceful)
    return Ray::INVALID_RAY_DATA_INDEX;

  const auto & other_map = aux ? _ray_data_map : _ray_aux_data_map;
  if (other_map.find(name) != other_map.end())
    mooseError(_error_prefix,
               ": Ray data with name ",
               name,
               " was not found.\n",
               "However, Ray ",
               (aux ? "(non-aux)" : "aux"),
               " data with said name was found.\n",
               "Did you mean to use ",
               (aux ? "getRayAuxDataIndex()/getRayAuxDataIndices()?"
                    : "getRayDataIndex()/getRayDataIndices()"));

  mooseError(_error_prefix, ": Unknown Ray ", (aux ? "aux" : ""), " data with name ", name);
}

std::vector<RayDataIndex>
RayTracingStudy::getRayDataIndicesInternal(const std::vector<std::string> & names,
                                           const bool aux,
                                           const bool graceful) const
{
  std::vector<RayDataIndex> indices(names.size());
  for (std::size_t i = 0; i < names.size(); ++i)
    indices[i] = getRayDataIndexInternal(names[i], aux, graceful);
  return indices;
}

const std::string &
RayTracingStudy::getRayDataNameInternal(const RayDataIndex index, const bool aux) const
{
  Threads::spin_mutex::scoped_lock lock(_spin_mutex);

  if (aux)
  {
    if (rayAuxDataSize() >= index)
      mooseError(_error_prefix, ": Unknown Ray aux data with index ", index);
  }
  else if (rayDataSize() >= index)
    mooseError(_error_prefix, ": Unknown Ray data with index ", index);

  return aux ? _ray_aux_data_names[index] : _ray_data_names[index];
}

std::vector<std::string>
RayTracingStudy::getRayDataNamesInternal(const std::vector<RayDataIndex> & indices,
                                         const bool aux) const
{
  std::vector<std::string> names(indices.size());
  for (std::size_t i = 0; i < indices.size(); ++i)
    names[i] = getRayDataNameInternal(indices[i], aux);
  return names;
}

RayDataIndex
RayTracingStudy::registerRayData(const std::string & name)
{
  return registerRayDataInternal(name, /* aux = */ false);
}

std::vector<RayDataIndex>
RayTracingStudy::registerRayData(const std::vector<std::string> & names)
{
  return registerRayDataInternal(names, /* aux = */ false);
}

RayDataIndex
RayTracingStudy::getRayDataIndex(const std::string & name, const bool graceful /* = false */) const
{
  return getRayDataIndexInternal(name, /* aux = */ false, graceful);
}

std::vector<RayDataIndex>
RayTracingStudy::getRayDataIndices(const std::vector<std::string> & names,
                                   const bool graceful /* = false */) const
{
  return getRayDataIndicesInternal(names, /* aux = */ false, graceful);
}

const std::string &
RayTracingStudy::getRayDataName(const RayDataIndex index) const
{
  return getRayDataNameInternal(index, /* aux = */ false);
}

std::vector<std::string>
RayTracingStudy::getRayDataNames(const std::vector<RayDataIndex> & indices) const
{
  return getRayDataNamesInternal(indices, /* aux = */ false);
}

RayDataIndex
RayTracingStudy::registerRayAuxData(const std::string & name)
{
  return registerRayDataInternal(name, /* aux = */ true);
}

std::vector<RayDataIndex>
RayTracingStudy::registerRayAuxData(const std::vector<std::string> & names)
{
  return registerRayDataInternal(names, /* aux = */ true);
}

RayDataIndex
RayTracingStudy::getRayAuxDataIndex(const std::string & name,
                                    const bool graceful /* = false */) const
{
  return getRayDataIndexInternal(name, /* aux = */ true, graceful);
}

std::vector<RayDataIndex>
RayTracingStudy::getRayAuxDataIndices(const std::vector<std::string> & names,
                                      const bool graceful /* = false */) const
{
  return getRayDataIndicesInternal(names, /* aux = */ true, graceful);
}

const std::string &
RayTracingStudy::getRayAuxDataName(const RayDataIndex index) const
{
  return getRayDataNameInternal(index, /* aux = */ true);
}

std::vector<std::string>
RayTracingStudy::getRayAuxDataNames(const std::vector<RayDataIndex> & indices) const
{
  return getRayDataNamesInternal(indices, /* aux = */ true);
}

bool
RayTracingStudy::hasRayKernels(const THREAD_ID tid)
{
  std::vector<RayKernelBase *> result;
  getRayKernels(result, tid);
  return result.size();
}

void
RayTracingStudy::getRayKernels(std::vector<RayKernelBase *> & result, SubdomainID id, THREAD_ID tid)
{
  // If the cache doesn't have any attributes yet, it means that we haven't set
  // the conditions yet. We do this so that it can be generated on the fly on first use.
  if (!_threaded_cache_ray_kernel[tid].numAttribs())
  {
    if (!_called_initial_setup)
      mooseError(_error_prefix, "Should not call getRayKernels() before initialSetup()");

    auto query = _fe_problem.theWarehouse()
                     .query()
                     .condition<AttribRayTracingStudy>(this)
                     .condition<AttribSystem>("RayKernel")
                     .condition<AttribThread>(tid);
    _threaded_cache_ray_kernel[tid] = query.clone();
  }

  _threaded_cache_ray_kernel[tid].queryInto(result, id);
}

void
RayTracingStudy::getRayKernels(std::vector<RayKernelBase *> & result,
                               SubdomainID id,
                               THREAD_ID tid,
                               RayID ray_id)
{
  // No Ray registration: no need to sift through objectsK
  if (!_use_ray_registration)
  {
    getRayKernels(result, id, tid);
  }
  // Has Ray registration: only pick the objects associated with ray_id
  else
  {
    // Get all of the kernels on this block
    std::vector<RayKernelBase *> rkbs;
    getRayKernels(rkbs, id, tid);

    // The RayTracingObjects associated with this ray
    const auto & ray_id_rtos = _threaded_ray_object_registration[tid][ray_id];

    // The result is the union of all of the kernels and the objects associated with this Ray
    result.clear();
    for (auto rkb : rkbs)
      if (ray_id_rtos.count(rkb))
        result.push_back(rkb);
  }
}

void
RayTracingStudy::getRayBCs(std::vector<RayBoundaryConditionBase *> & result,
                           BoundaryID id,
                           THREAD_ID tid)
{
  // If the cache doesn't have any attributes yet, it means that we haven't set
  // the conditions yet. We do this so that it can be generated on the fly on first use.
  if (!_threaded_cache_ray_bc[tid].numAttribs())
  {
    if (!_called_initial_setup)
      mooseError(_error_prefix, "Should not call getRayBCs() before initialSetup()");

    auto query = _fe_problem.theWarehouse()
                     .query()
                     .condition<AttribRayTracingStudy>(this)
                     .condition<AttribSystem>("RayBoundaryCondition")
                     .condition<AttribThread>(tid);
    _threaded_cache_ray_bc[tid] = query.clone();
  }

  _threaded_cache_ray_bc[tid].queryInto(result, std::make_tuple(id, false));
}

void
RayTracingStudy::getRayBCs(std::vector<RayBoundaryConditionBase *> & result,
                           const std::vector<TraceRayBndElement> & bnd_elems,
                           THREAD_ID tid,
                           RayID ray_id)
{
  // No Ray registration: no need to sift through objects
  if (!_use_ray_registration)
  {
    if (bnd_elems.size() == 1)
      getRayBCs(result, bnd_elems[0]._bnd_id, tid);
    else
    {
      std::vector<BoundaryID> bnd_ids(bnd_elems.size());
      for (MooseIndex(bnd_elems.size()) i = 0; i < bnd_elems.size(); ++i)
        bnd_ids[i] = bnd_elems[i]._bnd_id;
      getRayBCs(result, bnd_ids, tid);
    }
  }
  // Has Ray registration: only pick the objects associated with ray_id
  else
  {
    // Get all of the RayBCs on these boundaries
    std::vector<RayBoundaryConditionBase *> rbcs;
    if (bnd_elems.size() == 1)
      getRayBCs(rbcs, bnd_elems[0]._bnd_id, tid);
    else
    {
      std::vector<BoundaryID> bnd_ids(bnd_elems.size());
      for (MooseIndex(bnd_elems.size()) i = 0; i < bnd_elems.size(); ++i)
        bnd_ids[i] = bnd_elems[i]._bnd_id;
      getRayBCs(rbcs, bnd_ids, tid);
    }

    // The RayTracingObjects associated with this ray
    const auto & ray_id_rtos = _threaded_ray_object_registration[tid][ray_id];

    // The result is the union of all of the kernels and the objects associated with this Ray
    result.clear();
    for (auto rbc : rbcs)
      if (ray_id_rtos.count(rbc))
        result.push_back(rbc);
  }
}

std::vector<RayTracingObject *>
RayTracingStudy::getRayTracingObjects()
{
  std::vector<RayTracingObject *> result;
  _fe_problem.theWarehouse().query().condition<AttribRayTracingStudy>(this).queryInto(result);
  return result;
}

Real
RayTracingStudy::getBankedRayData(const RayID ray_id, const RayDataIndex index) const
{
  if (!_use_ray_registration)
    mooseError(_error_prefix, ": Should not use getBankedRayData() with Ray registration diabled");
  if (!_bank_rays_on_completion)
    mooseError(_error_prefix,
               ": Cannot use getBankedRayData() with bank_rays_on_completion = false");
  if (index >= rayDataSize())
    mooseError(_error_prefix,
               ": While getting a banked Ray data, the index ",
               index,
               " was supplied but the max possible index is ",
               rayDataSize() - 1);

  const Ray * ray = nullptr;
  for (const auto & possible_ray : _ray_bank)
    if (possible_ray->id() == ray_id)
      ray = possible_ray.get();

  const Real invalid_value = -std::numeric_limits<Real>::max();
  Real value = ray ? ray->data(index) : invalid_value;
  _communicator.max(value);

  if (value == invalid_value)
    mooseError(_error_prefix,
               ": While getting a banked Ray value, the Ray with id ",
               ray_id,
               " was not found in the Ray banks");

  return value;
}

Real
RayTracingStudy::getBankedRayAuxData(const RayID ray_id, const RayDataIndex index) const
{
  if (!_use_ray_registration)
    mooseError(_error_prefix,
               ": Should not use getBankedRayAuxData() with Ray registration diabled");
  if (!_bank_rays_on_completion)
    mooseError(_error_prefix,
               ": Cannot use getBankedRayData() with bank_rays_on_completion = false");
  if (index >= rayAuxDataSize())
    mooseError(_error_prefix,
               ": While getting a banked Ray aux data, the index ",
               index,
               " was supplied but the max possible index is ",
               rayAuxDataSize() - 1);

  const Ray * ray = nullptr;
  for (const auto & possible_ray : _ray_bank)
    if (possible_ray->id() == ray_id)
      ray = possible_ray.get();

  const Real invalid_value = std::numeric_limits<Real>::min();
  Real value = ray ? ray->auxData(index) : invalid_value;
  _communicator.max(value);

  if (value == invalid_value)
    mooseError(_error_prefix,
               ": While getting a banked Ray aux data, the Ray with id ",
               ray_id,
               " was not found in the Ray banks");

  return value;
}

RayID
RayTracingStudy::registerRay(const std::string & name, const bool graceful /* = false */)
{
  Threads::spin_mutex::scoped_lock lock(_spin_mutex);

  if (!_use_ray_registration)
    mooseError(_error_prefix, ": Cannot use registerRay() with Ray registration diabled");

  // This is parallel only for now. We could likely stagger the ID building like we do with
  // the unique IDs, but it would require a sync point which isn't there right now
  libmesh_parallel_only(comm());

  const auto search = _registered_ray_map.find(name);
  if (search != _registered_ray_map.end())
  {
    if (graceful)
      return Ray::INVALID_RAY_ID;
    else
      mooseError(_error_prefix,
                 ": Attempting to register Ray with name ",
                 name,
                 " but a Ray is already registered with that name!");
  }

  const auto next_id = _threaded_ray_object_registration[0].size();
  for (THREAD_ID tid = 0; tid < libMesh::n_threads(); ++tid)
    _threaded_ray_object_registration[tid].emplace_back();
  _registered_ray_map.emplace(name, next_id);
  _reverse_registered_ray_map.emplace(next_id, name);
  return next_id;
}

std::vector<RayID>
RayTracingStudy::registerRays(const std::vector<std::string> & names,
                              const bool graceful /* = false */)
{
  std::vector<RayID> ids(names.size());
  for (std::size_t i = 0; i < names.size(); ++i)
    ids[i] = registerRay(names[i], graceful);
  return ids;
}

RayID
RayTracingStudy::registeredRayID(const std::string & name, const bool graceful /* = false */) const
{
  Threads::spin_mutex::scoped_lock lock(_spin_mutex);

  if (!_use_ray_registration)
    mooseError(_error_prefix, ": Should not use registeredRayID() with Ray registration diabled");

  const auto search = _registered_ray_map.find(name);
  if (search != _registered_ray_map.end())
    return search->second;

  if (graceful)
    return Ray::INVALID_RAY_ID;

  mooseError(_error_prefix,
             ": Attempted to obtain ID of registered Ray ",
             name,
             ", but a Ray with said name is not registered.");
}

const std::string &
RayTracingStudy::registeredRayName(const RayID ray_id) const
{
  Threads::spin_mutex::scoped_lock lock(_spin_mutex);

  if (!_use_ray_registration)
    mooseError(_error_prefix, ": Should not use registeredRayName() with Ray registration diabled");

  const auto search = _reverse_registered_ray_map.find(ray_id);
  if (search != _reverse_registered_ray_map.end())
    return search->second;

  mooseError(_error_prefix,
             ": Attempted to obtain name of registered Ray with ID ",
             ray_id,
             ", but a Ray with said ID is not registered.");
}

Real
RayTracingStudy::computeTotalVolume()
{
  Real volume = 0;
  for (const auto & elem : *_mesh.getActiveLocalElementRange())
    volume += elem->volume();
  _communicator.sum(volume);
  return volume;
}

const std::vector<std::pair<unsigned short, std::vector<BoundaryID>>> *
RayTracingStudy::getInternalSidesets(const Elem * elem) const
{
  mooseAssert(_use_internal_sidesets, "Not using internal sidesets");

  if (!_has_internal_sidesets)
    return nullptr;

  mooseAssert(_mesh.getMesh().has_elem_integer(_internal_sideset_elem_integer_name),
              "Mesh doesn't have extra integer for indexing");

  // The offset for this element into _internal_sidesets;
  const auto offset = elem->get_extra_integer(_internal_sideset_elem_integer);

  // If it's invalid, it means this elem has no internal sidesets
  if (offset == DofObject::invalid_id)
    return nullptr;

  mooseAssert(_internal_sidesets.size() > offset, "Offset not in internal sidesets");
  return &_internal_sidesets[offset];
}

TraceData &
RayTracingStudy::initThreadedCachedTrace(const std::shared_ptr<Ray> & ray, THREAD_ID tid)
{
  mooseAssert(shouldCacheTrace(ray), "Not caching trace");
  mooseAssert(currentlyPropagating(), "Should only use while tracing");

  _threaded_cached_traces[tid]->emplace_back(ray);
  return _threaded_cached_traces[tid]->back();
}

void
RayTracingStudy::verifyRay(const std::shared_ptr<Ray> & ray,
                           const THREAD_ID tid,
                           const std::string & error_prefix)
{
  if (!ray)
    mooseError(error_prefix, ": Cannot verify null Ray");

  const auto error = [this, &error_prefix, &ray](const std::string & error) {
    mooseError(error_prefix, ":\n\nRay ", error, "\n\n", ray->getInfo(this));
  };

  if (ray->hasTraced())
    error("has uncleared counters");

  if (ray->invalidID())
    error("has an invalid ID");

  if (ray->invalidStartPoint())
    error("has an invalid start point");

  // This is unfortunately only a loose check for non-rectangular domains
  if (ray->endSet() &&
      !boundingBox().contains_point(ray->currentPoint() + ray->direction() * ray->maxDistance()))
    error("has an end point that is not within the mesh");

  if (!ray->startingElem() || !_mesh.queryElemPtr(ray->startingElem()->id()))
    error("has an unset starting element");

  if (!ray->startingElem()->active())
    error("has a starting element that is not active"
          "\n\nInactive starting elem information\n" +
          ray->startingElem()->get_info());

  if (!ray->startingElem()->contains_point(ray->startPoint()))
    error("has a starting element that does not contain the starting point"
          "\n\nStarting elem information\n" +
          ray->startingElem()->get_info());

  if (!ray->invalidStartingIncomingSide())
  {
    if (ray->startingElem()->n_sides() < ray->startingIncomingSide())
      error("has an starting incoming side that is not valid for its starting element"
            "\n\nStarting elem information\n" +
            ray->startingElem()->get_info());

    const Elem * side_elem = sidePtrHelper(ray->startingElem(), ray->startingIncomingSide(), tid);

    if (!side_elem->contains_point(ray->startPoint()))
      error("has an starting incoming side that does not contain the starting point"
            "\n\nStarting elem information\n" +
            ray->startingElem()->get_info() + "\n\nIncoming side information\n" +
            side_elem->get_info());

    if (!sideIsIncoming(ray->startingElem(), ray->startingIncomingSide(), ray->direction(), tid))
      error("has an incoming side that is not incoming\n\nCurrent elem information\n" +
            ray->startingElem()->get_info() + "\n\nIncoming side information\n" +
            side_elem->get_info());
  }

  if (ray->data().size() != rayDataSize())
    error("has data which is not sized as required by the study");
  if (ray->auxData().size() != rayAuxDataSize())
    error("has aux data which is not sized as required by the study");

  if (!ray->shouldContinue())
    error("should not continue");
}

void
RayTracingStudy::verifyUniqueRayIDs(const std::vector<std::shared_ptr<Ray>>::const_iterator begin,
                                    const std::vector<std::shared_ptr<Ray>>::const_iterator end,
                                    const bool global,
                                    const std::string & error_suffix) const
{
  auto print_id = [](const RayID id) {
    const std::string s = id == Ray::INVALID_RAY_ID ? "INVALID_RAY_ID" : std::to_string(id);
    return s;
  };

  // Check our local Rays first
  std::unordered_map<RayID, const Ray *> my_id_map;
  for (auto it = begin; it < end; ++it)
  {
    const std::shared_ptr<Ray> & ray = *it;
    if (!ray)
      continue;

    const auto emplace = my_id_map.emplace(ray->id(), ray.get());
    const bool found = !emplace.second;

    if (found)
    {
      const Ray * other_ray = emplace.first->second;
      mooseError(_error_prefix,
                 ":\nMultiple Rays exist with ID ",
                 print_id(ray->id()),
                 " on processor ",
                 _pid,
                 " ",
                 error_suffix,
                 "\n\nOffending Ray information:\n\n",
                 ray->getInfo(this),
                 "\n",
                 other_ray->getInfo(this));
    }
  }

  if (global)
  {
    // Package our IDs to send to all processors
    std::vector<RayID> my_ids;
    my_ids.reserve(my_id_map.size());
    for (const auto & id_ray_pair : my_id_map)
      my_ids.push_back(id_ray_pair.first);
    std::map<processor_id_type, std::vector<RayID>> send_ids;
    for (processor_id_type pid = 0; pid < n_processors(); ++pid)
      if (pid != _pid)
        send_ids.emplace(pid, my_ids);

    // Verify other processor's Rays against ours
    const auto check_ids = [this, &my_id_map, &error_suffix, &print_id](
                               processor_id_type pid, const std::vector<RayID> & ids) {
      for (const RayID id : ids)
      {
        const auto find = my_id_map.find(id);
        if (find != my_id_map.end())
        {
          const Ray * ray = find->second;
          mooseError(_error_prefix,
                     ":\nA Ray with ID ",
                     print_id(id),
                     " exists on processors ",
                     pid,
                     " and ",
                     _pid,
                     " ",
                     error_suffix,
                     "\n\nLocal offending Ray information:\n",
                     ray->getInfo(this));
        }
      }
    };
    Parallel::push_parallel_vector_data(_communicator, send_ids, check_ids);
  }
}

void
RayTracingStudy::verifyUniqueRays(const std::vector<std::shared_ptr<Ray>>::const_iterator begin,
                                  const std::vector<std::shared_ptr<Ray>>::const_iterator end,
                                  const std::string & error_suffix)
{
  std::set<const Ray *> rays;
  for (auto it = begin; it < end; ++it)
  {
    const std::shared_ptr<Ray> & ray = *it;

    if (!rays.insert(ray.get()).second) // false if not inserted into rays
      mooseError(_error_prefix,
                 ":\nMultiple shared_ptrs were found that point to the same Ray ",
                 error_suffix,
                 "\n\nOffending Ray:\n",
                 ray->getInfo(this));
  }
}

void
RayTracingStudy::moveRayToBuffer(std::shared_ptr<Ray> & ray,
                                 const THREAD_ID tid,
                                 std::string verify_ray_error_prefix)
{
  if (!ray)
    return;

#ifdef NDEBUG
  if (_verify_rays) // only guard with this in optimized modes
#endif
  {
    if (verify_ray_error_prefix == "")
      verify_ray_error_prefix = "While moving Ray to buffer in " + _error_prefix;
    verifyRay(ray, tid, verify_ray_error_prefix);
  }

  _parallel_ray_study->moveObjectToBuffer(ray, tid);
}

void
RayTracingStudy::moveRaysToBuffer(std::vector<std::shared_ptr<Ray>> & rays,
                                  const THREAD_ID tid,
                                  std::string verify_ray_error_prefix)
{
#ifdef NDEBUG
  if (_verify_rays) // only guard with this in optimized modes
#endif
  {
    if (verify_ray_error_prefix == "")
      verify_ray_error_prefix = "While moving Ray to buffer in " + _error_prefix;

    for (const std::shared_ptr<Ray> & ray : rays)
      if (ray)
        verifyRay(ray, tid, verify_ray_error_prefix);
  }

  _parallel_ray_study->moveObjectsToBuffer(rays, tid);
}

void
RayTracingStudy::reserveRayBuffer(const std::size_t size)
{
  if (!currentlyGenerating())
    mooseError(_error_prefix, ": Can only reserve in Ray buffer during generateRays()");

  _parallel_ray_study->reserveBuffer(size);
}

const Point &
RayTracingStudy::getSideNormal(const Elem * elem, unsigned short side, const THREAD_ID tid)
{
  std::unordered_map<std::pair<const Elem *, unsigned short>, Point> & cache =
      _threaded_cached_normals[tid];

  // See if we've already cached this side normal
  const auto elem_side_pair = std::make_pair(elem, side);
  const auto search = cache.find(elem_side_pair);

  // Haven't cached this side normal: compute it and then cache it
  if (search == cache.end())
  {
    _threaded_fe_face[tid]->reinit(elem, side);
    const auto & normal = _threaded_fe_face[tid]->get_normals()[0];
    cache.emplace(elem_side_pair, normal);
    return normal;
  }

  // Have cached this side normal: simply return it
  return search->second;
}

bool
RayTracingStudy::sameLevelActiveElems() const
{
  unsigned int min_level = std::numeric_limits<unsigned int>::max();
  unsigned int max_level = std::numeric_limits<unsigned int>::min();

  for (const auto & elem : *_mesh.getActiveLocalElementRange())
  {
    const auto level = elem->level();
    min_level = std::min(level, min_level);
    max_level = std::max(level, max_level);
  }

  _communicator.min(min_level);
  _communicator.max(max_level);

  return min_level == max_level;
}

Real
RayTracingStudy::subdomainHmax(const SubdomainID subdomain_id) const
{
  const auto find = _subdomain_hmax.find(subdomain_id);
  if (find == _subdomain_hmax.end())
    mooseError(_error_prefix, ": Subdomain ", subdomain_id, " not found in subdomain hmax map");
  return find->second;
}

bool
RayTracingStudy::isRectangularDomain() const
{
  Real bbox_volume = 1;
  for (unsigned int d = 0; d < _mesh.dimension(); ++d)
    bbox_volume *= std::abs(_b_box.max()(d) - _b_box.min()(d));

  return MooseUtils::absoluteFuzzyEqual(bbox_volume, totalVolume(), TOLERANCE);
}

void
RayTracingStudy::resetUniqueRayIDs()
{
  Threads::spin_mutex::scoped_lock lock(_spin_mutex);

  if (currentlyGenerating() || currentlyPropagating())
    mooseError("Cannot initialize unique Ray IDs during Ray generation or propagation");

  for (THREAD_ID tid = 0; tid < libMesh::n_threads(); ++tid)
    _threaded_next_ray_id[tid] = (RayID)_pid * (RayID)libMesh::n_threads() + (RayID)tid;
}

RayID
RayTracingStudy::generateUniqueRayID(const THREAD_ID tid)
{
  // Get the current ID to return
  const auto id = _threaded_next_ray_id[tid];

  // Advance so that the next call has the correct ID
  _threaded_next_ray_id[tid] += (RayID)n_processors() * (RayID)libMesh::n_threads();

  return id;
}

void
RayTracingStudy::resetRayCounters(const std::shared_ptr<Ray> & ray)
{
  if (currentlyPropagating())
    mooseError(_error_prefix, ": Cannot reset Ray counters during ray tracing");

  ray->resetCounters({});
}

bool
RayTracingStudy::sideIsIncoming(const Elem * const elem,
                                const unsigned short side,
                                const Point & direction,
                                const THREAD_ID tid)
{
  const auto & normal = getSideNormal(elem, side, tid);
  const auto dot = normal * direction;
  return dot < TraceRayTools::TRACE_TOLERANCE;
}
