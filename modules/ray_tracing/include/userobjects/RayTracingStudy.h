//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneralUserObject.h"

// Local Includes
#include "ParallelRayStudy.h"
#include "RayTracingAttributes.h"
#include "SidePtrHelper.h"
#include "TraceData.h"
#include "TraceRayBndElement.h"

// MOOSE Includes
#include "TheWarehouse.h"

// libMesh includes
#include "libmesh/mesh.h"

// Forward Declarations
class RayBoundaryConditionBase;
class RayKernelBase;
class TraceRay;
class RayTracingObject;

/**
 * Base class for Ray tracing studies that will generate Rays and then propoagate
 * all of them to termination.
 *
 * Subclasses _must_ override generateRays()
 */
class RayTracingStudy : public GeneralUserObject, public SidePtrHelper
{
public:
  RayTracingStudy(const InputParameters & parameters);

  static InputParameters validParams();

  virtual void initialSetup() override;
  virtual void residualSetup() override;
  virtual void jacobianSetup() override;
  virtual void subdomainSetup() override {}
  virtual void meshChanged() override;
  virtual void timestepSetup() override;

  /**
   * Setup for on subdomain change or subdomain AND ray change during ray tracing
   * @param subdomain The subdomain changed to
   * @param tid Thread id
   * @param ray_id ID of the ray initiating the change
   */
  virtual void subdomainSetup(SubdomainID subdomain, THREAD_ID tid, RayID ray_id);
  /**
   * Reinitialize objects for a Ray segment for ray tracing
   * @param elem The elem the segment is in
   * @param start Start point of the segment
   * @param end End point of the segment
   * @param length The length of the start -> end segment
   * @param tid Thread id
   */
  virtual void reinitSegment(
      const Elem * elem, const Point & start, const Point & end, const Real length, THREAD_ID tid);

  /**
   * Called at the end of a Ray segment
   * @param tid Thread id
   */
  virtual void postOnSegment(THREAD_ID tid);

  /**
   * Method for executing the study so that it can be called out of the standard UO execute()
   */
  void executeStudy();

  /**
   * Total number of processor crossings for Rays that finished on this processor
   */
  unsigned long long int endingProcessorCrossings() const { return _ending_processor_crossings; }
  /**
   * Max number of total processor crossings for Rays that finished on this processor
   */
  unsigned int endingMaxProcessorCrossings() const { return _ending_max_processor_crossings; }
  /**
   * Total number of processor crossings
   */
  unsigned long long int totalProcessorCrossings() const { return _total_processor_crossings; }
  /**
   * Max number of processor crossings for all Rays
   */
  unsigned int maxProcessorCrossings() const { return _max_processor_crossings; }

  /**
   * Total number of Ray/element intersections for Rays that finished on this processor
   */
  unsigned long long int endingIntersections() const { return _ending_intersections; }
  /**
   * Max number of intersections for Rays that finished on this processor
   */
  unsigned int endingMaxIntersections() const { return _ending_max_intersections; }
  /**
   * Total number of Ray/element intersections
   */
  unsigned long long int totalIntersections() const { return _total_intersections; }
  /**
   * Max number of intersections for a Ray
   */
  unsigned int maxIntersections() const { return _max_intersections; }
  /**
   * Max number of trajectory changes for a Ray
   */
  unsigned int maxTrajectoryChanges() const { return _max_trajectory_changes; }

  /**
   * Total amount of distance traveled by the rays that end on this processor
   */
  Real endingDistance() const { return _ending_distance; }
  /**
   * Total distance traveled by all Rays
   */
  Real totalDistance() const { return _total_distance; }

  unsigned long long int localTraceRayResult(const int result) const
  {
    return _local_trace_ray_results[result];
  }

  const ParallelRayStudy & parallelRayStudy() const { return *_parallel_ray_study; }

  /**
   * Max distance any Ray can travel
   */
  Real rayMaxDistance() const { return _ray_max_distance; }

  /**
   * Duration for execute() in seconds
   */
  Real executionTime() { return std::chrono::duration<Real>(_execution_time).count(); }
  /**
   * Duration for execute() in nanoseconds
   */
  Real executionTimeNano()
  {
    return std::chrono::duration<Real, std::nano>(_execution_time).count();
  }
  /**
   * Duration for creation of all Rays in seconds
   */
  Real generationTime() const { return std::chrono::duration<Real>(_generation_time).count(); }
  /**
   * Duration for creation of all Rays in seconds
   */
  Real propagationTime() const { return std::chrono::duration<Real>(_propagation_time).count(); }

  /**
   * Whether or not to tolerate failure
   */
  bool tolerateFailure() const { return _tolerate_failure; }

  /**
   * Whether or not to bank Rays on completion
   */
  bool bankRaysOnCompletion() const { return _bank_rays_on_completion; }

  /**
   * Whether or not to use Ray dependent subdomain setup
   */
  bool rayDependentSubdomainSetup() const { return _ray_dependent_subdomain_setup; }

  /**
   * Register a value to be filled in the data on a Ray with a given name.
   * @param name The name to register
   * \return The Ray data index for the registered value
   *
   * Note that this does not actually allocate the storage. It is simply a registration system to
   * keep track of the values stored in the Ray data.
   */
  RayDataIndex registerRayData(const std::string & name);
  /**
   * Register values to be filled in the data on a Ray with a given name.
   * @param names The names to register
   * \return The Ray data indices for the registered value
   *
   * Note that this does not actually allocate the storage. It is simply a registration system to
   * keep track of the values stored in the Ray data.
   */
  std::vector<RayDataIndex> registerRayData(const std::vector<std::string> & names);
  /**
   * Gets the index associated with a registered value in the Ray data
   * @param name The value name to get the index of
   * @param graceful Whether or not to exit gracefully if none is found (return
   * INVALID_RAY_DATA_INDEX)
   * \return The index for the value
   */
  RayDataIndex getRayDataIndex(const std::string & name, const bool graceful = false) const;
  /**
   * Gets the indices associated with registered values in the Ray data
   * @param names The value names to get the indices of
   * @param graceful Whether or not to exit gracefully if none is found (index is filled with
   * INVALID_RAY_DATA_INDEX)
   * \return The indices for the values
   */
  std::vector<RayDataIndex> getRayDataIndices(const std::vector<std::string> & names,
                                              const bool graceful = false) const;
  /**
   * Gets the name associated with a registered value in the Ray data
   * @param index The index to get the name of
   * \return The name associated with index
   */
  const std::string & getRayDataName(const RayDataIndex index) const;
  /**
   * Gets the names associated with registered values in the Ray data
   * @param indices The indices to get the names of
   * \return The associated names
   */
  std::vector<std::string> getRayDataNames(const std::vector<RayDataIndex> & indices) const;

  /**
   * The registered size of values in the Ray data.
   */
  std::size_t rayDataSize() const { return _ray_data_names.size(); }
  /**
   * Whether or not any Ray data are registered
   */
  bool hasRayData() const { return _ray_data_names.size(); }
  /**
   * The Ray data names
   */
  const std::vector<std::string> & rayDataNames() const { return _ray_data_names; }

  /**
   * Register a value to be filled in the aux data on a Ray with a given name.
   * @param name The name to register
   * \returns The Ray aux data index for the registered value
   *
   * Note that this does not actually allocate the storage. It is simply a registration system to
   * keep track of the values stored in the Ray aux data.
   */
  RayDataIndex registerRayAuxData(const std::string & name);
  /**
   * Register values to be filled in the aux data on a Ray with a given name.
   * @param names The names to register
   * \returns The Ray aux data indices for the registered values
   *
   * Note that this does not actually allocate the storage. It is simply a registration system to
   * keep track of the values stored in the Ray aux data.
   */
  std::vector<RayDataIndex> registerRayAuxData(const std::vector<std::string> & names);
  /**
   * Gets the index associated with a registered value in the Ray aux data
   * @param name The value name to get the index of
   * @param graceful Whether or not to exit gracefully if none is found (return
   * INVALID_RAY_DATA_INDEX)
   * \returns The index for the value
   */
  RayDataIndex getRayAuxDataIndex(const std::string & name, const bool graceful = false) const;
  /**
   * Gets the indices associated with registered values in the Ray aux data
   * @param names The value names to get the indices of
   * @param graceful Whether or not to exit gracefully if none is found (index is filled with
   * INVALID_RAY_DATA_INDEX)
   * \return The indices for the values
   */
  std::vector<RayDataIndex> getRayAuxDataIndices(const std::vector<std::string> & names,
                                                 const bool graceful = false) const;
  /**
   * Gets the name associated with a registered value in the Ray aux data
   * @param index The index to get the name of
   * \return The name associated with index
   */
  const std::string & getRayAuxDataName(const RayDataIndex index) const;
  /**
   * Gets the names associated with registered values in the Ray aux data
   * @param indices The indices to get the names of
   * \return The associated names
   */
  std::vector<std::string> getRayAuxDataNames(const std::vector<RayDataIndex> & indices) const;
  /**
   * The registered size of values in the Ray aux data.
   */
  std::size_t rayAuxDataSize() const { return _ray_aux_data_names.size(); }
  /**
   * Whether or not any Ray aux data are registered
   */
  bool hasRayAuxData() const { return _ray_aux_data_names.size(); }
  /**
   * The Ray aux data names
   */
  const std::vector<std::string> & rayAuxDataNames() const { return _ray_aux_data_names; }

  /// Hardware counters
  long long _cycles_count, _instructions_count, _clock_count, _cache_miss_count, _cache_load_count;

  /**
   * Whether or not there are currently any active RayKernel objects
   */
  bool hasRayKernels(const THREAD_ID tid);
  /**
   * Fills the active RayKernels associated with this study and a block into result
   */
  void getRayKernels(std::vector<RayKernelBase *> & result, SubdomainID id, THREAD_ID tid);
  /**
   * Fills the active RayKernels associated with this study into result
   */
  template <typename T>
  void getRayKernels(std::vector<T *> & result, THREAD_ID tid)
  {
    _fe_problem.theWarehouse()
        .query()
        .condition<AttribRayTracingStudy>(this)
        .condition<AttribSystem>("RayKernel")
        .condition<AttribThread>(tid)
        .queryInto(result);
  }
  /**
   * Fills th active RayKernels associeted with this study, block, and potentially Ray into result
   */
  void
  getRayKernels(std::vector<RayKernelBase *> & result, SubdomainID id, THREAD_ID tid, RayID ray_id);
  /**
   * Fills the active RayBCs associated with this study and a boundary into result
   */
  void getRayBCs(std::vector<RayBoundaryConditionBase *> & result, BoundaryID id, THREAD_ID tid);
  /**
   * Fills the active RayBCs associated with this study and boundaries result
   */
  template <typename T>
  void getRayBCs(std::vector<T *> & result, const std::vector<BoundaryID> & ids, THREAD_ID tid)
  {
    _fe_problem.theWarehouse()
        .query()
        .condition<AttribRayTracingStudy>(this)
        .condition<AttribSystem>("RayBoundaryCondition")
        .condition<AttribBoundaries>(ids)
        .condition<AttribThread>(tid)
        .queryInto(result);
  }
  /**
   * Fills the active RayBCs associated with this study into result
   */
  template <typename T>
  void getRayBCs(std::vector<T *> & result, THREAD_ID tid)
  {
    _fe_problem.theWarehouse()
        .query()
        .condition<AttribRayTracingStudy>(this)
        .condition<AttribSystem>("RayBoundaryCondition")
        .condition<AttribThread>(tid)
        .queryInto(result);
  }
  /**
   * Fills the active RayBCs associated with thie study, boundary elements, and potentially Ray into
   * result.
   *
   * This is purposely virtual because it allows derived studies to optimize the retrevial of RayBCs
   * during the trace in TraceRay.
   */
  virtual void getRayBCs(std::vector<RayBoundaryConditionBase *> & result,
                         const std::vector<TraceRayBndElement> & bnd_elems,
                         THREAD_ID tid,
                         RayID ray_id);

  /**
   * Gets all of the currently active RayTracingObjects
   */
  std::vector<RayTracingObject *> getRayTracingObjects();

  /**
   * Gets the current RayKernels for a thread, which are set in subdomainSetup()
   */
  const std::vector<RayKernelBase *> & currentRayKernels(THREAD_ID tid) const
  {
    return _threaded_current_ray_kernels[tid];
  }

  /**
   * Get the nodal bounding box for the domain
   */
  const BoundingBox & boundingBox() const { return _b_box; }
  /**
   * Get the inflated maximum length across the domain
   */
  Real domainMaxLength() const { return _domain_max_length; }
  /**
   * Get the current total volume of the domain
   */
  Real totalVolume() const { return _total_volume; }
  /**
   * Whether or not the domain is rectangular (if it is prefectly encompassed by its bounding box)
   */
  bool isRectangularDomain() const;

  /**
   * Whether or not the mesh has internal sidesets that ALSO have RayBCs on them
   *
   * NOTE: if useInternalSidesets() == false, this will be false even if the mesh HAS internal
   * sidesets
   */
  bool hasInternalSidesets() const { return _has_internal_sidesets; }
  /**
   * Get the internal sideset pairs (side : boundary ids) for internal sidesets that have an active
   * RayBC on them for an element
   * @param elem The element
   * @return Pointer to the pairs (nullptr if none)
   */
  const std::vector<std::pair<unsigned short, std::vector<BoundaryID>>> *
  getInternalSidesets(const Elem * elem) const;
  /**
   * Whether or not the mesh has active elements of the same level
   *
   * Use this over sameLevelActiveElems(), which is for internally setting
   * _has_same_level_active_elems
   */
  bool hasSameLevelActiveElems() const { return _has_same_level_active_elems; }

  /**
   * Acquire a Ray from the pool of Rays on thread \p tid.
   *
   * You should NEVER build Ray objects on your own. If you want a new Ray object,
   * you should always use this method. It allows for the re-use of Ray objects
   * after their destructors are called.
   *
   * This method is thread safe.
   *
   * \p args allows for the proper forwarding of arguments to the Ray constructor.
   * If reusing an object from the pool, the Ray::reset() method will be
   * called with said forwarded arguemnts.
   */
  template <typename... Args>
  std::shared_ptr<Ray> acquireRay(const THREAD_ID tid, Args &&... args)
  {
    return _parallel_ray_study->acquireObject(tid, std::forward<Args>(args)...);
  }

  /**
   * Access to the libMesh MeshBase.
   *
   * This is needed for unpack routines for a Ray, which has a context of this study.
   */
  MeshBase & meshBase() const { return _mesh; }

  /**
   * Get the outward normal for a given element side.
   */
  virtual const Point &
  getSideNormal(const Elem * elem, const unsigned short side, const THREAD_ID tid);
  /**
   * Gets the outward normals for a given element.
   * Returns a pointer to the normal for the zeroth side.
   */
  virtual const Point * getElemNormals(const Elem * /* elem */, const THREAD_ID /* tid */)
  {
    mooseError("Not implemented");
  }

  /**
   * Gets the data value for a banked ray with a given ID
   */
  Real getBankedRayData(const RayID ray_id, const RayDataIndex index) const;
  /**
   * Gets the data value for a banked ray with a given ID
   */
  Real getBankedRayAuxData(const RayID ray_id, const RayDataIndex index) const;

  /**
   * Gets the ID of a registered ray
   * @param name The name of said ray
   * @param graceful Whether or not to exit gracefully if none is found (with invalid_id)
   */
  RayID registeredRayID(const std::string & name, const bool graceful) const;
  /**
   * Gets the name of a registered ray
   * @param ray_id The ID of said ray
   */
  const std::string & registeredRayName(const RayID ray_id) const;

  /**
   * Whether or not ray registration is being used
   */
  bool useRayRegistration() const { return _use_ray_registration; }

  /**
   * Whether or not the residual is currently being computed
   */
  bool computingResidual() const { return _computing_residual; }
  /**
   * Whether or not the Jacobian is currently being computed
   */
  bool computingJacobian() const { return _computing_jacobian; }

  /**
   * Whether or not to store the Ray data on the cached Ray traces
   */
  bool dataOnCacheTraces() const { return _data_on_cache_traces; }
  /**
   * Whether or not to store the Ray aux data on the cached Ray traces
   */
  bool auxDataOnCacheTraces() const { return _aux_data_on_cache_traces; }
  /**
   *  Whether or not to cache individual element segments when _cache_traces = true
   */
  bool segmentsOnCacheTraces() const { return _segments_on_cache_traces; }

  /**
   * Virtual that allows for selection in if a Ray should be cached or not (only used when
   * _cache_traces).
   */
  virtual bool shouldCacheTrace(const std::shared_ptr<Ray> & /* ray */) const
  {
    return _always_cache_traces;
  }

  /**
   * Initialize a Ray in the threaded cached trace map to be filled with segments
   */
  TraceData & initThreadedCachedTrace(const std::shared_ptr<Ray> & ray, THREAD_ID tid);
  /**
   * Get the cached trace data structure
   */
  std::shared_ptr<std::vector<TraceData>> getCachedTraces() const { return _cached_traces; }

  /**
   * Get the cached hmax for all elements in a subdomain
   *
   * Used for scaling tolerances in ray tracing.
   */
  Real subdomainHmax(const SubdomainID subdomain_id) const;

  /**
   * Entry point for acting on a ray when it is completed (shouldContinue() == false)
   */
  virtual void onCompleteRay(const std::shared_ptr<Ray> & ray);

  /**
   * Adds a ray to the buffer to be traced.
   *
   * During generateRays(), this method can ONLY be called on thread 0. It is NOT thread
   * safe during generateRays()!
   *
   * During the propagation of Rays (accessing Rays within RayBC and RayKernel objects during
   * tracing), this method is thread safe.
   *
   * This moves the Ray into the buffer (with std::move), therefore \p ray will be nullptr after
   * this call.
   *
   * @param rays The Rays to move to the buffer
   * @param tid The thread id
   * @param verify_ray_error_prefix Error prefix to pass to verifyRay() to add context to the error.
   * If none is given, will default to a general error with study information
   */
  void moveRayToBuffer(std::shared_ptr<Ray> & ray,
                       const THREAD_ID tid = 0,
                       std::string verify_ray_error_prefix = "");
  /**
   * Adds rays to the buffer to be traced.
   *
   * During generateRays(), this method can ONLY be called on thread 0. It is NOT thread
   * safe during generateRays()!
   *
   * During the propagation of Rays (accessing Rays within RayBC and RayKernel objects during
   * tracing), this method is thread safe.
   *
   * This moves the rays into the buffer (with std::move), therefore all valid rays
   * in \p rays will be nullptr after this call.
   *
   * @param rays The Rays to move to the buffer
   * @param tid The thread id
   * @param verify_ray_error_prefix Error prefix to pass to verifyRay() to add context to the error.
   * If none is given, will default to a general error with study information
   */
  void moveRaysToBuffer(std::vector<std::shared_ptr<Ray>> & rays,
                        const THREAD_ID tid = 0,
                        std::string verify_ray_error_prefix = "");

  /**
   * Generates a unique RayID to be used for a Ray.
   *
   * Be sure if you are using this during tracing (in a RayKernel or a RayBC) to pass
   * in the correct thread id!
   *
   * By default, the unique IDs are never reset. If you wish to reset these such that
   * the IDs start at the beginning of the range again, say between calls to ray-tracing,
   * use resetUniqueRayIDs().
   */
  RayID generateUniqueRayID(const THREAD_ID tid = 0);

  /**
   * Verifies that the given Ray has the valid data needed to be traced.
   *
   * If used during generateRays(), it will also check the local uniqueness
   * of the Ray and the local uniqueness of the Ray's ID.
   *
   * @param ray The Ray to verify
   * @param tid The thread id
   * @param error_prefix Information to add at the beginning of the error if the
   * verification fails to add context
   */
  void verifyRay(const std::shared_ptr<Ray> & ray,
                 const THREAD_ID tid,
                 const std::string & error_prefix);

  /**
   * Verifies that the Rays in the given range have unique Ray IDs.
   *
   * @param begin The beginning of the range of Rays to check
   * @param end The end of the range of Rays to check
   * @param global If true, this will complete the verification across all processors
   * @param error_suffix Entry point for additional information in the error message
   */
  void verifyUniqueRayIDs(const std::vector<std::shared_ptr<Ray>>::const_iterator begin,
                          const std::vector<std::shared_ptr<Ray>>::const_iterator end,
                          const bool global,
                          const std::string & error_suffix) const;

  /**
   * Verifies that the Rays in the given range are unique. That is, that there are not multiple
   * shared_ptrs that point to the same Ray
   *
   * @param begin The beginning of the range of Rays to check
   * @param end The end of the range of Rays to check
   * @param error_suffix Entry point for additional information in the error message
   */
  void verifyUniqueRays(const std::vector<std::shared_ptr<Ray>>::const_iterator begin,
                        const std::vector<std::shared_ptr<Ray>>::const_iterator end,
                        const std::string & error_suffix);

  /**
   * Whether or not the study is propagating (tracing Rays)
   */
  bool currentlyPropagating() const { return _parallel_ray_study->currentlyExecuting(); }
  /**
   * Whether or not the study is generating
   */
  bool currentlyGenerating() const { return _parallel_ray_study->currentlyPreExecuting(); }

  /**
   * Casts the RayTracingStudy found in the given input parameters to a study of type T
   * with a meaningful error message if it fails
   *
   * Can be used for casting to derived studies on other objects that use them
   * (RayKernels, RayBCs, etc)
   */
  template <typename T>
  static T & castFromStudy(const InputParameters & params)
  {
    static_assert(std::is_base_of<RayTracingStudy, T>::value, "Not derived from a RayTracingStudy");

    RayTracingStudy * study =
        params.getCheckedPointerParam<RayTracingStudy *>("_ray_tracing_study");

    T * other_study = dynamic_cast<T *>(study);
    if (!other_study)
      ::mooseError(params.get<std::string>("_type"),
                   " '",
                   params.get<std::string>("_object_name"),
                   "' must be paired with a ",
                   typeid(T).name());
    return *other_study;
  }

  /**
   * Gets the threaded TraceRay object for \p tid.
   */
  const TraceRay & traceRay(const THREAD_ID tid) const { return *_threaded_trace_ray[tid]; }

  /**
   * Whether or not Ray verification is enabled in opt mode
   */
  bool verifyRays() const { return _verify_rays; }

  /**
   * Whether or not \p side is incoming on element \p elem in direction \p direction.
   */
  bool sideIsIncoming(const Elem * const elem,
                      const unsigned short side,
                      const Point & direction,
                      const THREAD_ID tid);

protected:
  /**
   * Subclasses should override this to determine how to generate Rays.
   * This will be called within execute() and makes up the "generation phase"
   * of the algorithm.
   */
  virtual void generateRays() = 0;

  /**
   * Entry point before study execution
   */
  virtual void preExecuteStudy();
  /**
   * Entry point after study execution
   */
  virtual void postExecuteStudy();

  /**
   * Registers a Ray with a given name.
   * @param name The name to register
   * @param graceful Whether or not to error on registration failure. If true, Ray::INVALID_RAY_ID
   * will be returned.
   * @return The allocated ID for the registered ray
   */
  RayID registerRay(const std::string & name, const bool graceful = false);
  /**
   * Registers Rays with the given names
   * @param names The names to register
   * @param graceful Whether or not to error on registration failure. If true, Ray::INVALID_RAY_ID
   * will be returned.
   * @return The allocated IDs for said registered rays
   */
  std::vector<RayID> registerRays(const std::vector<std::string> & name,
                                  const bool graceful = false);

  /**
   * Helper function for computing the total domain volume
   */
  Real computeTotalVolume();

  /**
   * Gets the threaded TraceRay objects.
   *
   * Allows for other studies to change options in the TraceRay objects.
   */
  const std::vector<std::shared_ptr<TraceRay>> & threadedTraceRay() const
  {
    return _threaded_trace_ray;
  }

  /**
   * Gets the writeable current RayKernels for a thread
   *
   * Allows for other ray studies to fill the current ray kernels in a custom manner
   */
  std::vector<RayKernelBase *> & currentRayKernelsWrite(THREAD_ID tid)
  {
    return _threaded_current_ray_kernels[tid];
  }

  /**
   * Reserve \p size entires in the Ray buffer.
   *
   * This can only be used within generateRays() and should be used when possible
   * to reserve entires in the buffer before adding Rays via moveRay(s)ToBuffer.
   */
  void reserveRayBuffer(const std::size_t size);

  /**
   * Determine whether or not the mesh currently has active elements that are
   * all the same level
   */
  bool sameLevelActiveElems() const;

  /**
   * Builds quadrature points for a given segment using the _segment_qrule
   * @param start Start point of the segment
   * @param end End point of the segment
   * @param length The lengh of the start -> end segment
   * @param points Points to fill into (should be sized ahead of time)
   * @param weights Weights to fill into (should be sized ahead of time)
   */
  virtual void buildSegmentQuadrature(const Point & start,
                                      const Point & end,
                                      const Real length,
                                      std::vector<Point> & points,
                                      std::vector<Real> & weights) const;

#ifndef NDEBUG
  /**
   * Get the debug ray bank.
   *
   * This is only available outside of optimized mode. All Rays that finish on this processor
   * are inserted into this bank for debugging purposes.
   */
  const std::vector<std::shared_ptr<Ray>> & debugRayBank() const { return _debug_ray_bank; }
#endif

  /**
   * Resets a Ray's counters.
   *
   * This is the only way to reset a Ray's counters other than obtaining a Ray from the pool
   * or constructing a new one. This is done on purpose in order to disallow the resetting
   * of any of a Ray's counters during ray-tracing, which results in undefined behavior
   * in the trace routines.
   */
  void resetRayCounters(const std::shared_ptr<Ray> & ray);

  /**
   * Resets the generation of unique RayIDs via generateUniqueRayID() to the beginning
   * of the range.
   */
  void resetUniqueRayIDs();

  /// The Mesh
  MooseMesh & _mesh;
  /// The Communicator
  const Parallel::Communicator & _comm;
  /// The rank of this processor (this actually takes time to lookup - so just do it once)
  const processor_id_type _pid;
  /// Prefix used in errors (type() 'name()':)
  const std::string _error_prefix;

  /// Whether or not to perform coverage checks on RayKernels
  const bool _ray_kernel_coverage_check;
  /// Whether not to check if all element faces are planar
  const bool _planar_face_check;
  /// Whether or not to use Ray registration
  const bool _use_ray_registration;
  /// Whether or not to use the internal sidesets in ray tracing
  const bool _use_internal_sidesets;
  /// Whether or not to tolerate a Ray Tracing failure
  const bool _tolerate_failure;
  /// Whether or not to bank rays on completion
  const bool _bank_rays_on_completion;
  /// Whether or not subdomain setup is dependent on the Ray
  const bool _ray_dependent_subdomain_setup;

  /// Whether or not to cache traces on every trace execution
  bool _always_cache_traces;
  /// Whether or not to store the Ray data on the cache traces
  const bool _data_on_cache_traces;
  /// Whether or not to store the Ray aux data on the cache traces
  const bool _aux_data_on_cache_traces;
  /// Whether or not to cache individual element segments when caching
  const bool _segments_on_cache_traces;
  /// Max distance a Ray can travel before being killed (can change)
  const Real _ray_max_distance;
  /// Whether or not to verify Rays that are added to the buffer in opt mode
  const bool _verify_rays;

private:
  /**
   * Perform coverage checks (coverage of RayMaterials and RayKernels, if enabled)
   */
  void coverageChecks();

  /**
   * Perform checks to see if the listed dependencies in the RayTracingObjects exist
   */
  void dependencyChecks();
  /**
   * Verifies that the dependencies exist for a set of RayTracingObjects
   */
  void verifyDependenciesExist(const std::vector<RayTracingObject *> & rtos);

  /**
   * Check for if all of the element types in the mesh are supported by ray tracing
   */
  void traceableMeshChecks();

  /**
   * Does the setup for internal sidesets. This includes:
   * - Setting if the mesh has internal sidesets (_has_internal_sidesets)
   * - Setting up the _internal_sidesets data structure for getInternalSidesets() getter
   *
   * Note: if _use_internal_sidesets = false, even if the mesh has internal sidesets, no setup
   * will be done and _has_internal_sidesets will be false.
   */
  void internalSidesetSetup();

  /**
   * Sets up the maps from Ray to associated RayTracingObjects if _use_ray_registration
   */
  void associateRegisteredRays();

  /**
   * Zero the AuxVariables that the registered AuxRayKernels contribute to.
   */
  void zeroAuxVariables();

  /**
   * Caches the hmax for all elements in each subdomain
   */
  void setupSubdomainHmax();

  /**
   * These don't mean anything here:
   */
  virtual void execute() override;
  virtual void initialize() override {}
  virtual void finalize() override {}

  /**
   * Internal method for registering Ray data or Ray aux data with a name.
   */
  RayDataIndex registerRayDataInternal(const std::string & name, const bool aux);
  /**
   * Internal method for registering Ray data or Ray aux data with names.
   */
  std::vector<RayDataIndex> registerRayDataInternal(const std::vector<std::string> & names,
                                                    const bool aux);
  /**
   * Internal method for getting the index of Ray data or Ray aux data.
   */
  RayDataIndex
  getRayDataIndexInternal(const std::string & name, const bool aux, const bool graceful) const;
  /**
   * Internal method for getting the indicies of Ray data or Ray aux data.
   */
  std::vector<RayDataIndex> getRayDataIndicesInternal(const std::vector<std::string> & names,
                                                      const bool aux,
                                                      const bool graceful) const;

  /**
   * Internal method for getting the name of Ray data or Ray aux data.
   */
  const std::string & getRayDataNameInternal(const RayDataIndex index, const bool aux) const;
  /**
   * Internal method for getting the names of Ray data or Day aux data.
   */
  std::vector<std::string> getRayDataNamesInternal(const std::vector<RayDataIndex> & indices,
                                                   const bool aux) const;

  /// Timing
  //@{
  PerfID _execute_study_timer;
  PerfID _generate_timer;
  PerfID _propagate_timer;
  std::chrono::steady_clock::time_point _execution_start_time;
  std::chrono::steady_clock::duration _execution_time;
  std::chrono::steady_clock::duration _generation_time;
  std::chrono::steady_clock::duration _propagation_time;
  //@}

  /// The map from Ray data names to index
  std::unordered_map<std::string, RayDataIndex> _ray_data_map;
  /// The map from Ray aux data names to index
  std::unordered_map<std::string, RayDataIndex> _ray_aux_data_map;

  /// The names for each Ray data entry
  std::vector<std::string> _ray_data_names;
  /// The names for each Ray aux data entry
  std::vector<std::string> _ray_aux_data_names;

  /// Map from registered Ray name to ID
  std::unordered_map<std::string, RayID> _registered_ray_map;
  /// Map from registered Ray ID to name
  std::unordered_map<RayID, std::string> _reverse_registered_ray_map;

  /// Storage for the cached traces
  std::shared_ptr<std::vector<TraceData>> _cached_traces;
  /// The threaded storage for cached traces
  std::vector<std::unique_ptr<std::vector<TraceData>>> _threaded_cached_traces;

  /// Whether or not the Jacobian is currently being computed (set in jacobianSetup())
  bool _computing_jacobian;
  /// Whether or not the residual is currently being computed (set in residualSetup())
  bool _computing_residual;
  /// Number of currently cached objects for Jacobian/residual for each thread
  std::vector<unsigned int> _num_cached;

  /// Whether or not the mesh has internal sidesets
  bool _has_internal_sidesets;
  /// The name to use for the extra element integer that indexes into _internal_sidesets
  const std::string _internal_sideset_elem_integer_name;
  /// The index into the extra elem integers for the value that indexes into _internal_sidesets
  unsigned int _internal_sideset_elem_integer;
  /// Internal sideset data, if any (indexed with the extra elem integer named _internal_sideset_elem_integer_name)
  std::vector<std::vector<std::pair<unsigned short, std::vector<BoundaryID>>>> _internal_sidesets;
  /// Whether or not the mesh has active elements of the same level
  bool _has_same_level_active_elems;

  /// Nodal bounding box for the domain
  BoundingBox _b_box;
  /// An inflated max distance for the domain
  Real _domain_max_length;
  /// The total volume of the domain
  Real _total_volume;

  /// Threaded cached subdomain query for RayKernelBase objects pretaining to this study
  std::vector<TheWarehouse::QueryCache<AttribSubdomains>> _threaded_cache_ray_kernel;
  /// Threaded cached boundary query for RayBC objects pretaining to this study
  std::vector<TheWarehouse::QueryCache<AttribBoundaries>> _threaded_cache_ray_bc;
  /// Threaded storage for all of the RayTracingObjects associated with a single Ray
  std::vector<std::vector<std::set<const RayTracingObject *>>> _threaded_ray_object_registration;
  /// The current RayKernel objects for each thread
  std::vector<std::vector<RayKernelBase *>> _threaded_current_ray_kernels;
  /// The TraceRay objects for each thread (they do the physical tracing)
  std::vector<std::shared_ptr<TraceRay>> _threaded_trace_ray;
  /// Face FE used for computing face normals for each thread
  std::vector<std::unique_ptr<FEBase>> _threaded_fe_face;
  /// Face quadrature used for computing face normals for each thread
  std::vector<std::unique_ptr<QBase>> _threaded_q_face;
  /// Threaded cache for side normals that have been computed already during tracing
  std::vector<std::unordered_map<std::pair<const Elem *, unsigned short>, Point>>
      _threaded_cached_normals;
  /// Cumulative Ray bank
  std::vector<std::shared_ptr<Ray>> _ray_bank;
#ifndef NDEBUG
  /// All of the Rays that have ended on this processor for use in debugging
  std::vector<std::shared_ptr<Ray>> _debug_ray_bank;
#endif
  /// Storage for the next available unique RayID, obtained via generateUniqueRayID()
  std::vector<RayID> _threaded_next_ray_id;

  /// The study that used to actually execute (trace) the Rays
  const std::unique_ptr<ParallelRayStudy> _parallel_ray_study;

  /// Quadrature rule for laying points across a 1D ray segment
  UniquePtr<QBase> _segment_qrule;

  /// Total number of processor crossings for Rays that finished on this processor
  unsigned long long int _ending_processor_crossings;
  /// Max number of total processor crossings for Rays that finished on this processor
  unsigned int _ending_max_processor_crossings;
  /// Total number of processor crossings
  unsigned long long int _total_processor_crossings;
  /// Max number of processor crossings for all Rays
  unsigned int _max_processor_crossings;

  /// Total number of Ray/element intersections for Rays that finished on this processor
  unsigned long long int _ending_intersections;
  /// Max number of intersections for Rays that finished on this processor
  unsigned int _ending_max_intersections;
  /// Max number of trajectory changes for Rays that finished on this processor
  unsigned int _ending_max_trajectory_changes;
  /// Total number of Ray/element intersections
  unsigned long long int _total_intersections;
  /// Max number of intersections for single Ray
  unsigned int _max_intersections;
  /// Max number of trajectory changes for a single Ray
  unsigned int _max_trajectory_changes;

  /// Total distance traveled by Rays that end on this processor
  Real _ending_distance;
  /// Total distance traveled by all Rays
  Real _total_distance;

  /// Cumulative results on this processor from the threaded TraceRay objects
  std::vector<unsigned long long int> _local_trace_ray_results;

  /// The cached hmax for all elements in a subdomain
  std::unordered_map<SubdomainID, Real> _subdomain_hmax;

  /// Whether or not we've called initial setup - used to stop from late registration
  bool _called_initial_setup;

  /// Spin mutex object for locks
  mutable Threads::spin_mutex _spin_mutex;
};
