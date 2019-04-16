//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef PlanarAverage_H
#define PlanarAverage_H

#include "SpatialAverageBase.h"

class PlanarAverage;

template <>
InputParameters validParams<PlanarAverage>();

/**
 * Compute a planar average, of a variable as a function of distance
 * from the origin projected onto the axis vector. These averages
 * are effectively taken over planes orthogonal to axis.
 * A special case are planar averages along x, y, z axes.
 */
class PlanarAverage : public SpatialAverageBase
{
public:
  PlanarAverage(const InputParameters & parameters);

protected:
  /// compute the distance of the current quadarature point for binning
  virtual Real computeDistance() override;

  /// vector used for measuring distance
  const Point _axis;

  /// norm of the _axis vector
  const Real _axis_norm;
};

#endif // PlanarAverage_H
