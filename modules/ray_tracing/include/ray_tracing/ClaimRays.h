//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// Local includes
#include "Ray.h"
#include "SidePtrHelper.h"

// System includes
#include <unordered_map>

// libMesh includes
#include "libmesh/bounding_box.h"
#include "libmesh/point_locator_base.h"

// Forward declarations
class RayTracingStudy;
class MooseMesh;

/**
 * Helper object for claiming Rays
 */
class ClaimRays : public SidePtrHelper
{
public:
  /**
   * Constructor.
   * @param study The RayTracingStudy
   * @param mesh The MooseMesh
   * @param rays The vector of Rays that need to be claimed
   * @param local_rays Insertion point for Rays that have been claimed
   * @param do_exchange Whether or not an exhange is needed, i.e., if "rays" still needs to be
   * filled by objects on other processors
   */
  ClaimRays(RayTracingStudy & study,
            MooseMesh & mesh,
            const std::vector<std::shared_ptr<Ray>> & rays,
            std::vector<std::shared_ptr<Ray>> & local_rays,
            const bool do_exchange);

  /**
   * Initialize the object
   */
  void init();
  /**
   * Claim the Rays
   */
  void claim();

protected:
  /**
   * Entry point before claim()
   */
  virtual void preClaim() {}
  /**
   * Entry point after claim()
   */
  virtual void postClaim() {}
  /**
   * Entry point before possibly claiming a Ray
   */
  virtual void prePossiblyClaim(const std::shared_ptr<Ray> & /* ray */) {}
  /**
   * Entry point for acting on a Ray after it is claimed
   */
  virtual void postClaim(std::shared_ptr<Ray> & ray, const Elem * elem);

  /// The mesh
  MooseMesh & _mesh;
  /// The communicator
  const libMesh::Parallel::Communicator & _comm;
  /// This processor ID
  const processor_id_type _pid;

  /// Whether or not the Rays need to be initially exchanged
  const bool _do_exchange;

  /// The study, used for receive context in communicating a Ray
  RayTracingStudy & _study;

private:
  /**
   * Builds the bounding boxes (_bbox, _inflated_bboxes, _inflated_neighbor_bboxes).
   */
  void buildBoundingBoxes();

  /**
   * Build the map of elements to all of their point neighbors
   *
   * TODO: Move this eventually into MooseMesh, MeshBase, or FEProblemBase
   */
  void buildPointNeighbors();

  /**
   * Possibly claim a Ray.
   */
  void possiblyClaim(const std::shared_ptr<Ray> & obj);

  /**
   * Try to claim a spatial point.
   *
   * @param point The point to claim
   * @param id An ID associated with the point
   * @param elem The local element to first consider for this processor's ownership
   * @return The element that contains the point if we claim the point, nullptr if we don't claim it
   */
  const Elem * claimPoint(const Point & point, const RayID id, const Elem * elem);

  /**
   * Creates an inflated bounding box.
   */
  BoundingBox inflateBoundingBox(const BoundingBox & bbox, const Real multiplier = 0.01);

  /// The Rays that need to be searched to possibly claimed
  const std::vector<std::shared_ptr<Ray>> & _rays;
  /// The local Rays that are claimed
  std::vector<std::shared_ptr<Ray>> & _local_rays;

  /// The point locator
  std::unique_ptr<PointLocatorBase> _point_locator = nullptr;

  /// The bounding box for this processor
  BoundingBox _bbox;
  /// The global bounding box
  BoundingBox _global_bbox;
  /// The inflaed bounding boxes for all processors
  std::vector<BoundingBox> _inflated_bboxes;
  /// Inflated bounding boxes that are neighboring to this processor (pid : bbox for each entry)
  std::vector<std::pair<processor_id_type, BoundingBox>> _inflated_neighbor_bboxes;

  /// Map of point neighbors for each element
  std::unordered_map<dof_id_type, std::vector<const Elem *>> _elem_point_neighbors;
};
