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

/**
 * MeshGenerator for
 */
class LayeredElementIDGenerator : public MeshGenerator
{
public:
  static InputParameters validParams();

  LayeredElementIDGenerator(const InputParameters & parameters);

  std::unique_ptr<MeshBase> generate() override;

protected:
  /// helper to find layer given coordinate value
  unsigned int layerHelper(Real x) const;

  /// mesh to add the subdomain to
  std::unique_ptr<MeshBase> & _input;

  /// the direction normal to the layers
  unsigned int _direction;

  /// the layer bounds
  std::vector<Real> _bounds;

  /// the number of layers
  unsigned int _n_layer;
};
