//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PlanarAverage.h"

#include "libmesh/quadrature.h"

registerMooseObject("MooseApp", PlanarAverage);

template <>
InputParameters
validParams<PlanarAverage>()
{
  InputParameters params = validParams<SpatialAverageBase>();
  params.addRequiredParam<Point>("axis", "Vector used to measure distance from origin.");
  params.addClassDescription(
      "Compute a planar average, of a variable as a function of distance from the origin projected "
      "onto the axis vector. These averages are effectively taken over planes orthogonal to axis. "
      "A special case are planar averages along x, y, z axes.");
  return params;
}

PlanarAverage::PlanarAverage(const InputParameters & parameters)
  : SpatialAverageBase(parameters), _axis(getParam<Point>("axis")), _axis_norm(_axis.norm())
{
}

Real
PlanarAverage::computeDistance()
{
  return (_q_point[_qp] - _origin) * _axis / _axis_norm;
}
