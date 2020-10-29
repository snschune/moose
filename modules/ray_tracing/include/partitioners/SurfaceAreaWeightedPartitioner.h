//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE includes
#include "MooseEnum.h"
#include "PetscExternalPartitioner.h"

/**
 * Partitions a mesh with weights determined by element surface areas.
 */
class SurfaceAreaWeightedPartitioner : public PetscExternalPartitioner
{
public:
  SurfaceAreaWeightedPartitioner(const InputParameters & params);

  static InputParameters validParams();

  virtual std::unique_ptr<Partitioner> clone() const override;
  virtual dof_id_type computeElementWeight(Elem & elem) override;
  virtual dof_id_type computeSideWeight(Elem & elem, unsigned int side) override;

protected:
  virtual void _do_partition(MeshBase & mesh, const unsigned int n) override;

  Real _min_side_area;
  Real _min_elem_surface_area;
};
