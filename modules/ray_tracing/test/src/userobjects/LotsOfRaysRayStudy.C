//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LotsOfRaysRayStudy.h"

// Local includes
#include "RayTracingAngularQuadrature.h"

// libMesh includes
#include "libmesh/parallel_algebra.h"

registerMooseObject("RayTracingTestApp", LotsOfRaysRayStudy);

InputParameters
LotsOfRaysRayStudy::validParams()
{
  auto params = RepeatableRayStudyBase::validParams();

  params.addParam<bool>(
      "vertex_to_vertex",
      true,
      "Enable generation of rays from vertices on boundary sides in the direction of "
      "the vertices on the other side of the same elem");
  params.addParam<bool>("centroid_to_vertex",
                        true,
                        "Enable generation of rays from centroids on boundary sides in the "
                        "direction of all other nodes in the same elem");
  params.addParam<bool>("centroid_to_centroid",
                        true,
                        "Enable generation of rays from centroids on boundary sides to centroids "
                        "of all boundary elements");
  params.addParam<bool>("edge_to_edge",
                        false,
                        "Enable generation of rays from centroids of boundary edges in the "
                        "direction of all other edge centroids in the same elem");
  params.addParam<bool>("side_aq",
                        false,
                        "Enable generation of rays from boundary side centroids in the direction "
                        "of angular quadrature");
  params.addParam<bool>("centroid_aq",
                        false,
                        "Enable generation of rays from boundary element centroids in the "
                        "direction of angular quadrature.");

  params.addRangeCheckedParam<unsigned int>(
      "polar_quad_order",
      5,
      "polar_quad_order % 2",
      "Order of the polar quadrature angular quadrature generation. Polar angle is between ray and "
      "the normal direction. Must be odd.");
  params.addRangeCheckedParam<unsigned int>(
      "azimuthal_quad_order",
      2,
      "azimuthal_quad_order > 0",
      "Order of the azimuthal quadrature per quadrant for angular quadrature generation. The "
      "azimuthal angle is measured in a plane perpendicular to the normal. Not needed in 2D.");

  params.addParam<bool>(
      "compute_expected_distance",
      false,
      "Whether or not to compute the expected distance for all of the Rays generated");

  // When we define rays, they are already localized on the processors that start them
  // and the starting elem/incoming side is set
  params.set<bool>("_claim_after_define_rays") = false;
  // Can't have replicated Rays if the Rays are on their starting processors only
  params.set<bool>("_define_rays_replicated") = false;
  // Don't need to use registration here. For the purposes of testing and this study,
  // we'll really only ever access the Rays by ID (if at all)
  params.set<bool>("_use_ray_registration") = false;

  return params;
}

LotsOfRaysRayStudy::LotsOfRaysRayStudy(const InputParameters & parameters)
  : RepeatableRayStudyBase(parameters),
    _vertex_to_vertex(getParam<bool>("vertex_to_vertex")),
    _centroid_to_vertex(getParam<bool>("centroid_to_vertex")),
    _centroid_to_centroid(getParam<bool>("centroid_to_centroid")),
    _edge_to_edge(getParam<bool>("edge_to_edge")),
    _side_aq(getParam<bool>("side_aq")),
    _centroid_aq(getParam<bool>("centroid_aq")),
    _compute_expected_distance(getParam<bool>("compute_expected_distance")),
    _polar_quad_order(getParam<unsigned int>("polar_quad_order")),
    _azimuthal_quad_order(getParam<unsigned int>("azimuthal_quad_order"))
{
  // Computing the expected distance only works with rectanglar domains because we have to compute
  // the intersection with the Ray and the domain bounding box (which only represents the domain if
  // it is rectangular)
  if (_compute_expected_distance && !isRectangularDomain())
    mooseError(_error_prefix, " does not support non-rectangular domains");

  const bool polar_set = parameters.isParamSetByUser("polar_quad_order");
  const bool azimuthal_set = parameters.isParamSetByUser("azimuthal_quad_order");
  if (_mesh.dimension() == 1)
  {
    if (polar_set)
      paramError("polar_quad_order", "Not needed for 1D");
    if (azimuthal_set)
      paramError("azimuthal_quad_order", "Not needed for 1D");
  }
  if (!_side_aq && !_centroid_aq)
  {
    if (polar_set)
      paramError("polar_quad_order", "Not needed without side_aq or centroid_aq enabled");
    if (azimuthal_set)
      paramError("azimuthal_quad_order", "Not needed without side_aq or centroid_aq enabled");
  }
  if (azimuthal_set && _mesh.dimension() == 2)
    paramError("azimuthal_quad_order", "Not needed for 2D");
}

void
LotsOfRaysRayStudy::defineRays()
{
  if (_compute_expected_distance)
  {
    _bbox_intersection_helper =
        libmesh_make_unique<BoundingBoxIntersectionHelper>(boundingBox(), _mesh.dimension());
    _expected_distance = 0;
  }

  // Gather the boundary elem centroids
  std::vector<Point> boundary_centroids;
  if (_centroid_to_centroid)
  {
    for (const auto & elem : *_mesh.getActiveLocalElementRange())
      if (elem->on_boundary())
        boundary_centroids.push_back(elem->centroid());
    _comm.allgather(boundary_centroids);
  }

  std::vector<std::pair<Real, Real>> half_aq_angles;
  if (_side_aq || _centroid_aq)
  {
    std::vector<Real> unused_weights;
    if (_mesh.dimension() == 2)
      RayTracingAngularQuadrature::getHalfRange2D(
          _polar_quad_order, half_aq_angles, unused_weights, 0.5 * M_PI);
    else
      RayTracingAngularQuadrature::getHalfRange3D(
          4 * _azimuthal_quad_order, _polar_quad_order, half_aq_angles, unused_weights, 0.5 * M_PI);
  }

  for (const auto & bnd_elem : *_mesh.getBoundaryElementRange())
  {
    const auto side = bnd_elem->_side;
    const Elem * elem = bnd_elem->_elem;

    // Starting from local elems only
    if (elem->processor_id() != _pid)
      continue;

    // Centroid on this boundary face
    const auto & side_centroid = sidePtrHelper(elem, side)->centroid();

    // Vertices on this boundary side -> in direction of all other vertices not on said boundary
    // side on elem
    if (_vertex_to_vertex)
      for (const auto n : elem->nodes_on_side(side))
        if (elem->is_vertex(n))
          for (unsigned int v_to = 0; v_to < elem->n_vertices(); ++v_to)
            if (!elem->is_node_on_side(v_to, side))
              defineRay(elem, side, elem->point(n), elem->point(v_to), false);

    // Centroid on this boundary side -> in direction of all other vertices not on said boundary
    // side on elem
    if (_centroid_to_vertex)
      for (unsigned int v_to = 0; v_to < elem->n_vertices(); ++v_to)
        if (!elem->is_node_on_side(v_to, side))
          defineRay(elem, side, side_centroid, elem->point(v_to), false);

    // Centroid on this boundary side -> centroids of all boundary elements
    if (_centroid_to_centroid)
      for (const auto & other_centroid : boundary_centroids)
        defineRay(elem, RayTracingCommon::invalid_side, side_centroid, other_centroid, true);

    // Centroids of edges on this boundary side -> in direction of all other edge centroids
    if (_edge_to_edge)
      for (const auto edge : elem->edge_index_range())
        if (elem->is_edge_on_side(edge, side))
        {
          const Point edge_centroid = elem->build_edge_ptr(edge)->centroid();
          for (const auto edge_to : elem->edge_index_range())
            if (edge != edge_to)
              defineRay(
                  elem, side, edge_centroid, elem->build_edge_ptr(edge_to)->centroid(), false);
        }

    if (_side_aq)
    {
      const auto inward_normal = -1.0 * getSideNormal(elem, side, /* tid = */ 0);
      const auto rotation_matrix =
          RayTracingAngularQuadrature::getRotationMatrix(inward_normal, _mesh.dimension());

      for (std::size_t l = 0; l < half_aq_angles.size(); ++l)
        defineRay(elem,
                  side,
                  side_centroid,
                  side_centroid + RayTracingAngularQuadrature::getDirection(
                                      l, rotation_matrix, half_aq_angles, _mesh.dimension()),
                  false);
    }
  }

  if (_centroid_aq)
  {
    const auto dummy_rotation_matrix =
        RayTracingAngularQuadrature::getRotationMatrix(Point(1, 0, 0), _mesh.dimension());
    for (const Elem * elem : *_mesh.getActiveLocalElementRange())
      if (elem->processor_id() == _pid && elem->on_boundary())
      {
        const auto centroid = elem->centroid();

        for (std::size_t l = 0; l < half_aq_angles.size(); ++l)
        {
          const auto direction = RayTracingAngularQuadrature::getDirection(
              l, dummy_rotation_matrix, half_aq_angles, _mesh.dimension());

          defineRay(elem, RayTracingCommon::invalid_side, centroid, centroid + direction, false);
          defineRay(elem, RayTracingCommon::invalid_side, centroid, centroid - direction, false);
        }
      }
  }

  if (_compute_expected_distance)
    _communicator.sum(_expected_distance);

  // Insertion point for other dervied test objects to modify the Rays
  modifyRays();
}

void
LotsOfRaysRayStudy::defineRay(const Elem * starting_elem,
                              const unsigned short incoming_side,
                              const Point & p1,
                              const Point & p2,
                              const bool ends_within_mesh)
{
  std::shared_ptr<Ray> ray = acquireRay(/* tid = */ 0);

  ray->setStartingElem(starting_elem);
  ray->setStartingIncomingSide(incoming_side);
  ray->setID(generateUniqueRayID());

  if (ends_within_mesh)
  {
    ray->setStartEnd(p1, p2);

    if (_compute_expected_distance)
      _expected_distance += (p2 - p1).norm();
  }
  else
  {
    ray->setStartDirection(p1, p2 - p1);

    if (_compute_expected_distance)
    {
      const Point direction = (p2 - p1) / (p2 - p1).norm();
      const Point end = _bbox_intersection_helper->intersection(p1, direction);
      if (end == RayTracingCommon::invalid_point)
        mooseError(_error_prefix,
                   ": Expected distance end intersection not found\n\n",
                   ray->getInfo(this));

      _expected_distance += (end - p1).norm();
    }
  }

  _rays.emplace_back(std::move(ray));
}

void
LotsOfRaysRayStudy::modifyRays()
{
}
