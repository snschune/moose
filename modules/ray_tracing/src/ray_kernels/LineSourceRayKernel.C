//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LineSourceRayKernel.h"

// MOOSE includes
#include "Function.h"

registerMooseObject("RayTracingApp", LineSourceRayKernel);

InputParameters
LineSourceRayKernel::validParams()
{
  InputParameters params = RayKernel::validParams();

  params.addClassDescription(
      "Demonstrates the multiple ways that scalar values can be introduced "
      "into RayKernels, e.g. (controllable) constants, functions, "
      "postprocessors, and data on Rays. Implements the weak form $(\\psi_i, -f)$.");

  params.addParam<Real>("value", 1.0, "Coefficient to multiply by the line source term");
  params.addParam<FunctionName>("function", "1", "A function that describes the line source");
  params.addParam<PostprocessorName>(
      "postprocessor", 1, "A postprocessor whose value is multiplied by the line source");
  params.addParam<std::vector<std::string>>(
      "ray_data_factor_names", "The names of the Ray data to scale the residual by (if any)");
  params.addParam<std::vector<std::string>>(
      "ray_aux_data_factor_names",
      "The names of the Ray aux data to scale the residual by (if any)");

  params.declareControllable("value");

  return params;
}

LineSourceRayKernel::LineSourceRayKernel(const InputParameters & params)
  : RayKernel(params),
    _scale(getParam<Real>("value")),
    _function(getFunction("function")),
    _postprocessor(getPostprocessorValue("postprocessor")),
    _ray_data_factor_indices(
        _study.getRayDataIndices(getParam<std::vector<std::string>>("ray_data_factor_names"))),
    _ray_aux_data_factor_indices(_study.getRayAuxDataIndices(
        getParam<std::vector<std::string>>("ray_aux_data_factor_names")))
{
}

Real
LineSourceRayKernel::computeQpResidual()
{
  Real factor = _scale * _postprocessor * _function.value(_t, _q_point[_qp]);

  // Scale by any Ray data if given
  if (_ray_data_factor_indices.empty())
    for (const auto index : _ray_data_factor_indices)
      factor *= currentRay()->data(index);
  // Scale by any Ray aux data if given
  if (_ray_aux_data_factor_indices.empty())
    for (const auto index : _ray_aux_data_factor_indices)
      factor *= currentRay()->auxData(index);

  return _test[_i][_qp] * -factor;
}
