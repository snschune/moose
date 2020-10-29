//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SurfaceAreaWeightedPartitioner.h"

#include "libmesh/elem.h"

registerMooseObject("RayTracingApp", SurfaceAreaWeightedPartitioner);

InputParameters
SurfaceAreaWeightedPartitioner::validParams()
{
  return PetscExternalPartitioner::validParams();
}

SurfaceAreaWeightedPartitioner::SurfaceAreaWeightedPartitioner(const InputParameters & params)
  : PetscExternalPartitioner(params)
{
}

void
SurfaceAreaWeightedPartitioner::_do_partition(MeshBase & mesh, const unsigned int n)
{
  // Loop over all of the elements and find the _smallest_ side area
  _min_side_area = std::numeric_limits<Real>::max();
  _min_elem_surface_area = std::numeric_limits<Real>::max();

  for (const auto & elem : mesh.active_local_element_ptr_range())
  {
    Real elem_surface_area = 0;

    for (const auto s : elem->side_index_range())
    {
      const auto side_area = elem->side_ptr(s)->volume();

      _min_side_area = std::min(side_area, _min_side_area);

      elem_surface_area += side_area;
    }

    _min_elem_surface_area = std::min(elem_surface_area, _min_elem_surface_area);
  }

  // Find the min over all procs
  _communicator.min(_min_side_area);
  _communicator.min(_min_elem_surface_area);

  // Call the base class
  PetscExternalPartitioner::_do_partition(mesh, n);
}

std::unique_ptr<Partitioner>
SurfaceAreaWeightedPartitioner::clone() const
{
  return libmesh_make_unique<SurfaceAreaWeightedPartitioner>(_pars);
}

dof_id_type
SurfaceAreaWeightedPartitioner::computeElementWeight(Elem & elem)
{
  Real elem_surface_area = 0;

  for (const auto s : elem.side_index_range())
    elem_surface_area += elem.side_ptr(s)->volume();

  return 100 * (elem_surface_area / _min_elem_surface_area);
}

dof_id_type
SurfaceAreaWeightedPartitioner::computeSideWeight(Elem & elem, unsigned int side)
{
  const auto side_elem = elem.side_ptr(side);

  return 100 * (side_elem->volume() / _min_side_area);
}
