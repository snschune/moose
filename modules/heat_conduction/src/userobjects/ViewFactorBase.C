//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ViewFactorBase.h"
#include "libmesh/quadrature.h"

#include <limits>

template <>
InputParameters
validParams<ViewFactorBase>()
{
  InputParameters params = validParams<SideUserObject>();
  params.addParam<Real>("view_factor_tol",
                        std::numeric_limits<Real>::max(),
                        "Tolerance for checking view factors. Default is to allow everything.");
  params.addParam<bool>("normalize_view_factor",
                        true,
                        "Determines if view factors are normalized to sum to one (consistent with "
                        "their definition).");
  params.addClassDescription(
      "A base class for automatic computation of view factors between sidesets.");
  return params;
}

ViewFactorBase::ViewFactorBase(const InputParameters & parameters)
  : SideUserObject(parameters),
    _n_sides(boundaryIDs().size()),
    _areas(_n_sides),
    _view_factor_tol(getParam<Real>("view_factor_tol")),
    _normalize_view_factor(getParam<bool>("normalize_view_factor"))
{
  // sizing the view factor array
  _view_factors.resize(_n_sides);
  for (auto & v : _view_factors)
    v.resize(_n_sides);

  // set up the map from the side id to the local index & side name to local index
  std::vector<BoundaryName> boundary_names = getParam<std::vector<BoundaryName>>("boundary");
  for (unsigned int j = 0; j < boundary_names.size(); ++j)
  {
    _side_id_index[_mesh.getBoundaryID(boundary_names[j])] = j;
    _side_name_index[boundary_names[j]] = j;
  }
}

Real
ViewFactorBase::getViewFactor(BoundaryID from_id, BoundaryID to_id) const
{
  auto from = _side_id_index.find(from_id);
  auto to = _side_id_index.find(to_id);
  if (from == _side_id_index.end())
    mooseError("Boundary id ",
               from_id,
               " with name ",
               _mesh.getBoundaryName(from_id),
               " not listed in boundary parameter.");

  if (to == _side_id_index.end())
    mooseError("Boundary id ",
               to_id,
               " with name ",
               _mesh.getBoundaryName(to_id),
               " not listed in boundary parameter.");

  return _view_factors[from->second][to->second];
}

Real
ViewFactorBase::getViewFactor(BoundaryName from_name, BoundaryName to_name) const
{
  auto from = _side_name_index.find(from_name);
  auto to = _side_name_index.find(to_name);
  if (from == _side_name_index.end())
    mooseError("Boundary id ",
               _mesh.getBoundaryID(from_name),
               " with name ",
               from_name,
               " not listed in boundary parameter.");

  if (to == _side_name_index.end())
    mooseError("Boundary id ",
               _mesh.getBoundaryID(to_name),
               " with name ",
               to_name,
               " not listed in boundary parameter.");

  return _view_factors[from->second][to->second];
}

void
ViewFactorBase::finalize()
{
  // do some communication before finalizing view_factors
  for (unsigned int i = 0; i < _n_sides; ++i)
    for (unsigned int j = 0; j < _n_sides; ++j)
      gatherSum(_view_factors[i][j]);

  finalizeViewFactor();
  checkAndNormalizeViewFactor();
}

void
ViewFactorBase::threadJoin(const UserObject & y)
{
  const ViewFactorBase & pps = static_cast<const UnobstructedPlanarViewFactor &>(y);
  for (unsigned int i = 0; i < _n_sides; ++i)
  {
    for (unsigned int j = 0; j < _n_sides; ++j)
      _view_factors[i][j] += pps._view_factors[i][j];
  }
  threadJoinViewFactor(y);
}

void
ViewFactorBase::checkAndNormalizeViewFactor()
{
  for (unsigned int from = 0; from < _n_sides; ++from)
  {
    Real s = 0;
    for (unsigned int to = 0; to < _n_sides; ++to)
      s += _view_factors[from][to];

    if (std::abs(1 - s) > _view_factor_tol)
      mooseError("View factor from boundary ", boundaryNames()[from], " add to ", s);

    if (_normalize_view_factor)
      for (unsigned int to = 0; to < _n_sides; ++to)
        _view_factors[from][to] /= s;
  }
}
