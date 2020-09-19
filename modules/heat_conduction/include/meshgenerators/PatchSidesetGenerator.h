//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MeshGenerator.h"
#include "libmesh/replicated_mesh.h"

// Forward declarations

/**
 * Subdivides a sidesets into smaller patches each of which is going
 * to be a new patch
 */
class PatchSidesetGenerator : public MeshGenerator
{
public:
  static InputParameters validParams();

  PatchSidesetGenerator(const InputParameters & parameters);

  std::unique_ptr<MeshBase> generate() override;

  class Bubble
  {
  public:
    // default constructor
    Bubble();

    // constructor initializing with a single seed elem
    Bubble(const Elem * seed);

    // recomputes the frontier (non-incremental)
    void computeFrontier(const std::set<const Elem *> & assigned_elems);

    // removes an element from the frontier
    void removeFromFrontier(const Elem * elem);

    // removes elem from frontier and adds it as member elem
    void addToMember(const Elem * elem, const std::set<const Elem *> & assigned_elems);

    // returns workload
    Real workload() const { return _workload; }

    // how many frontier elements are present
    unsigned int frontierSize() const { return _frontier.size(); }

    // returns the frontier element closest to the seed
    const Elem * closestFrontierElem() const;

    // access to frontier
    std::set<const Elem *> frontier() const { return _frontier; }

    // access to members
    std::set<const Elem *> members() const { return _member_elems; }

    // find bubble centroid point and return the element whose centroid is closest
    const Elem * bubbleCentroidElem() const;

  protected:
    // the workload is amount of work assigned to this bubble
    Real _workload;

    // the seed/current center element
    const Elem * _seed;

    // set of elements that are assigned to this bubble
    std::set<const Elem *> _member_elems;

    // set of elements that are adjacent to at least one member_elem
    std::set<const Elem *> _frontier;
  };

protected:
  /// returns the name of the _n_patches subdivisions derived from _sideset
  std::vector<BoundaryName> sidesetNameHelper(const std::string & base_name) const;

  /// find the element that is farthest away from the seeds
  const Elem * mostDistantElement(const std::vector<const Elem *> & seeds) const;

  /// grows bubbles from seeds
  void growBubbles(const std::vector<const Elem *> & seeds, std::vector<Bubble> & bubbles) const;

  Elem * boundaryElementHelper(MeshBase & mesh, libMesh::ElemType type) const;

  /// a function for implementing custom partitioning
  void partition(MeshBase & mesh) const;

  std::unique_ptr<MeshBase> & _input;

  /// dimensionality of the sidesets to partition
  unsigned int _dim;

  /// the number of patches that this sideset generator divides _sideset into
  unsigned int _n_patches;

  /// The sideset that will be subdivided
  const boundary_id_type & _sideset;

  /// the name of the partitioner being used
  MooseEnum _partitioner_name;

  /// number of elements of the boundary mesh
  dof_id_type _n_boundary_mesh_elems;
};
