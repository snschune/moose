//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LayeredElementIDGenerator.h"
#include "Conversion.h"
#include "MooseMeshUtils.h"
#include "CastUniquePointer.h"

#include "libmesh/elem.h"

#include <algorithm>

registerMooseObject("MooseApp", LayeredElementIDGenerator);

InputParameters
LayeredElementIDGenerator::validParams()
{
  InputParameters params = MeshGenerator::validParams();

  params.addRequiredParam<MeshGeneratorName>("input", "The mesh we want to modify");
  MooseEnum directions("x=0 y=1 z=2");

  params.addRequiredParam<MooseEnum>("direction", directions, "The direction of the layers.");
  params.addRequiredParam<std::vector<Real>>("bounds",
                                             "The 'bounding' positions of the layers i.e.: '0, "
                                             "1.2, 3.7, 4.2' will mean 3 layers between those "
                                             "positions.");
  params.addRequiredParam<ExtraElementIDName>("extra_id_name", "The name of the extra element id");
  params.addClassDescription("");
  return params;
}

LayeredElementIDGenerator::LayeredElementIDGenerator(const InputParameters & parameters)
  : MeshGenerator(parameters),
    _input(getMesh("input")),
    _direction(getParam<MooseEnum>("direction")),
    _bounds(getParam<std::vector<Real>>("bounds")),
    _n_layer(_bounds.size() - 1)
{
  // make sure bounds has at least two entries
  if (_bounds.size() < 2)
    paramError("bounds", "Must have at least two entries");

  // make sure bounds is ascending
  for (unsigned int j = 1; j < _bounds.size(); ++j)
    if (_bounds[j] <= _bounds[j - 1])
      paramError("bounds", "Must be ascending");
}

std::unique_ptr<MeshBase>
LayeredElementIDGenerator::generate()
{
  std::unique_ptr<MeshBase> mesh = std::move(_input);

  // add the extra element id
  unsigned int extra_id = mesh->add_elem_integer(getParam<ExtraElementIDName>("extra_id_name"));

  // Loop over the elements & get subdomains
  std::vector<std::set<SubdomainID>> layer_subdomain_ids(_n_layer);
  for (const auto & elem : mesh->element_ptr_range())
  {
    Real coord = elem->centroid()(_direction);
    Real layer_id = layerHelper(coord);
    layer_subdomain_ids[layer_id].insert(elem->subdomain_id());
  }

  // if this mesh is a parallel mesh, communication is necessary so everyone knows
  // all the subdomain ids
  for (unsigned int j = 0; j < _n_layer; ++j)
    comm().set_union(layer_subdomain_ids[j]);

  // create a mapping from (layer_id, subdomain_id) => extra element id
  unsigned int current_extra_id = 0;
  std::vector<std::map<SubdomainID, unsigned int>> extra_id_map(_n_layer);
  for (unsigned int j = 0; j < _n_layer; ++j)
  {
    // get set for a single layer
    const auto & layer_set = layer_subdomain_ids[j];

    // sort the entries in layer set
    std::vector<SubdomainID> layer_vector(layer_set.begin(), layer_set.end());
    std::sort(layer_vector.begin(), layer_vector.end());

    // store assignment in map
    for (unsigned int i = 0; i < layer_vector.size(); ++i)
    {
      extra_id_map[j][layer_vector[i]] = current_extra_id;
      ++current_extra_id;
    }
  }

  // assign extra element ids
  for (const auto & elem : mesh->element_ptr_range())
  {
    Real coord = elem->centroid()(_direction);
    Real layer_id = layerHelper(coord);
    SubdomainID sub_id = elem->subdomain_id();
    auto it = extra_id_map[layer_id].find(sub_id);
    if (it == extra_id_map[layer_id].end())
      mooseError("Subdomain ", sub_id, " does not exist in layer ", layer_id);
    elem->set_extra_integer(extra_id, it->second);
  }

  return dynamic_pointer_cast<MeshBase>(mesh);
}

unsigned int
LayeredElementIDGenerator::layerHelper(Real x) const
{
  if (x < _bounds[0])
  {
    mooseDoOnce(mooseWarning("Element centroid out of bounds"));
    return 0;
  }

  if (x > _bounds.back())
  {
    mooseDoOnce(mooseWarning("Element centroid out of bounds"));
    return _n_layer - 1;
  }

  for (unsigned int j = 1; j < _bounds.size(); ++j)
    if (x <= _bounds[j])
      return j - 1;

  mooseError("Should never get here.");
}
