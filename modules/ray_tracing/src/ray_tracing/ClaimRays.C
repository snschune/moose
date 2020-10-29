//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ClaimRays.h"

// Local includes
#include "RayTracingStudy.h"

// libMesh includes
#include "libmesh/elem.h"
#include "libmesh/parallel_algebra.h"
#include "libmesh/parallel_sync.h"
#include "libmesh/mesh_tools.h"

ClaimRays::ClaimRays(RayTracingStudy & study,
                     MooseMesh & mesh,
                     const std::vector<std::shared_ptr<Ray>> & rays,
                     std::vector<std::shared_ptr<Ray>> & local_rays,
                     const bool do_exchange)
  : SidePtrHelper(),
    _mesh(mesh),
    _comm(_mesh.getMesh().comm()),
    _pid(_comm.rank()),
    _do_exchange(do_exchange),
    _study(study),
    _rays(rays),
    _local_rays(local_rays)
{
}

void
ClaimRays::claim()
{
  preClaim();

  _local_rays.clear();

  // Grab the point locator
  _point_locator = _mesh.getMesh().sub_point_locator();
  _point_locator->enable_out_of_mesh_mode();

  // Exchange: filter Rays into processors that _may_ claim them
  std::unordered_map<processor_id_type, std::vector<std::shared_ptr<Ray>>> rays_to_send;
  if (_do_exchange)
    for (auto & ray : _rays)
      for (processor_id_type pid = 0; pid < _comm.size(); ++pid)
        if (_inflated_bboxes[pid].contains_point(ray->currentPoint()))
          rays_to_send[pid].push_back(ray);

  // Functor for possibly claiming a vector of Rays
  std::function<void(processor_id_type, const std::vector<std::shared_ptr<Ray>> &)> claim_functor =
      [&](processor_id_type /* pid */, const std::vector<std::shared_ptr<Ray>> & rays) {
        for (auto & ray : rays)
          possiblyClaim(ray);
      };

  // Send the relevant Rays to everyone and then claim
  if (_do_exchange)
    Parallel::push_parallel_packed_range(_comm, rays_to_send, &_study, claim_functor);
  // Already have the relevant Rays, just claim
  else
    claim_functor(_pid, _rays);

  postClaim();
}

void
ClaimRays::possiblyClaim(const std::shared_ptr<Ray> & ray)
{
  prePossiblyClaim(ray);

  const auto elem =
      claimPoint(ray->currentPoint(), ray->id(), (*_point_locator)(ray->currentPoint()));
  if (elem)
  {
    _local_rays.push_back(ray);
    postClaim(_local_rays.back(), elem);
  }
}

const Elem *
ClaimRays::claimPoint(const Point & point, const RayID id, const Elem * elem)
{
  if (elem)
  {
    // Looking for smallest (even ID Ray) or largest (odd ID Ray) elem id
    const bool smallest = id % 2 == 0;

    // Start with the element we found, as it is a valid candidate
    const Elem * extremum_elem = elem;

    // All point neighbors for this element
    mooseAssert(_elem_point_neighbors.count(elem->id()), "Not in point neighbor map");
    const auto & neighbors = _elem_point_neighbors.at(elem->id());

    // Find element that matches the extremum criteria
    for (const auto & neighbor : neighbors)
    {
      mooseAssert(neighbor->active(), "Inactive neighbor");

      if ((smallest && neighbor->id() < extremum_elem->id()) || // satisfies
          (!smallest && neighbor->id() > extremum_elem->id()))  // ...one of the id checks
        if (neighbor->contains_point(point))                    // and also contains the point
          extremum_elem = neighbor;
    }

    // Claim the object if we own the extremum elem
    if (extremum_elem->processor_id() == _pid)
    {
      mooseAssert(extremum_elem->active(), "Inactive element");
      return extremum_elem;
    }
  }

  return nullptr;
}

void
ClaimRays::postClaim(std::shared_ptr<Ray> & ray, const Elem * elem)
{
  mooseAssert(_mesh.queryElemPtr(elem->id()) == elem, "Mesh doesn't contain elem");
  mooseAssert(elem->active(), "Inactive element");

  // Set the starting element
  ray->setStartingElem(elem);

  // If the user set a side, see if it is valid (contains start and is entrant). If it is not valid,
  // invalidate everything associated with it
  auto invalidated_side = RayTracingCommon::invalid_side;
  if (!ray->invalidStartingIncomingSide() &&
      (!sidePtrHelper(elem, ray->startingIncomingSide())->contains_point(ray->startPoint()) ||
       !_study.sideIsIncoming(elem, ray->startingIncomingSide(), ray->direction(), 0)))
  {
    invalidated_side = ray->startingIncomingSide();
    ray->invalidateStartingIncomingSide();
  }

  // Incoming side is invalid, so see if we can find one (contains start and is entrant). If the
  // user originally set one and it is wrong, don't re-check that side
  if (ray->invalidStartingIncomingSide())
    for (const auto s : elem->side_index_range())
      if (s != invalidated_side && sidePtrHelper(elem, s)->contains_point(ray->startPoint()) &&
          _study.sideIsIncoming(elem, s, ray->direction(), 0))
      {
        ray->setStartingIncomingSide(s);
        break;
      }
}

void
ClaimRays::init()
{
  buildBoundingBoxes();
  buildPointNeighbors();
}

void
ClaimRays::buildBoundingBoxes()
{
  // Local bounding box
  _bbox = MeshTools::create_local_bounding_box(_mesh.getMesh());
  _global_bbox = _bbox;

  // Gather the bounding boxes of all processors
  std::vector<std::pair<Point, Point>> bb_points = {static_cast<std::pair<Point, Point>>(_bbox)};
  _comm.allgather(bb_points, true);
  _inflated_bboxes.resize(_comm.size());
  for (processor_id_type pid = 0; pid < _comm.size(); ++pid)
  {
    const BoundingBox pid_bbox = static_cast<BoundingBox>(bb_points[pid]);
    _inflated_bboxes[pid] = inflateBoundingBox(pid_bbox);
    _global_bbox.union_with(pid_bbox);
  }

  // Find intersecting (neighbor) bounding boxes
  _inflated_neighbor_bboxes.clear();
  for (processor_id_type pid = 0; pid < _comm.size(); ++pid)
  {
    // Skip this processor
    if (pid == _pid)
      continue;
    // Insert if the searched processor's bbox intersects my bbox
    const auto & pid_bbox = _inflated_bboxes[pid];
    if (_bbox.intersects(pid_bbox))
      _inflated_neighbor_bboxes.emplace_back(pid, pid_bbox);
  }
}

void
ClaimRays::buildPointNeighbors()
{
  _elem_point_neighbors.clear();
  const auto & node_to_elem_map = _mesh.nodeToElemMap();

  for (const auto & elem : _mesh.getMesh().active_element_ptr_range())
  {
    auto & fill = _elem_point_neighbors[elem->id()];
    for (unsigned int v = 0; v < elem->n_vertices(); ++v)
    {
      const auto & node = elem->node_ptr(v);
      for (const auto & neighbor_id : node_to_elem_map.at(node->id()))
      {
        if (neighbor_id == elem->id())
          continue;

        const auto & neighbor = _mesh.elemPtr(neighbor_id);
        if (std::count(fill.begin(), fill.end(), neighbor) == 0)
          fill.emplace_back(neighbor);
      }
    }
  }
}

BoundingBox
ClaimRays::inflateBoundingBox(const BoundingBox & bbox, const Real multiplier)
{
  Real amount = multiplier * (bbox.max() - bbox.min()).norm();
  Point inflation(amount, amount, amount);
  auto inflated_bbox = bbox;
  inflated_bbox.first -= inflation;
  inflated_bbox.second += inflation;
  return inflated_bbox;
}
