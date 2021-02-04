//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CoupledConvectiveLineSource.h"
#include "LayeredAverage.h"
// MOOSE includes
#include "Function.h"

registerMooseObject("RayTracingApp", CoupledConvectiveLineSource);
registerMooseObject("RayTracingApp", ADCoupledConvectiveLineSource);

template <bool is_ad>
InputParameters
CoupledConvectiveLineSourceTempl<is_ad>::validParams()
{
  InputParameters params = GenericRayKernel<is_ad>::validParams();

  params.addClassDescription("Blabla");

  params.addRequiredCoupledVar("fluid_temperature", "Temperature of the fluid...");
  params.addRequiredCoupledVar("htc", "Heat transfer coefficient...");
  params.addRequiredParam<FunctionName>("heated_perimeter", "The heated perimeter");

  return params;
}

template <bool is_ad>
CoupledConvectiveLineSourceTempl<is_ad>::CoupledConvectiveLineSourceTempl(const InputParameters & params)
  : GenericRayKernel<is_ad>(params),
    _Tf(this->template coupledGenericValue<is_ad>("fluid_temperature")),
    _htc(this->template coupledGenericValue<is_ad>("htc")),
    _heated_perimeter(getFunction("heated_perimeter"))
{
}

template <bool is_ad>
GenericReal<is_ad>
CoupledConvectiveLineSourceTempl<is_ad>::computeQpResidual()
{
  return _test[_i][_qp] * _htc[_qp] * _heated_perimeter.value(_t, _q_point[_qp]) * (_u[_qp] - _Tf[_qp]);
}

template class CoupledConvectiveLineSourceTempl<false>;
template class CoupledConvectiveLineSourceTempl<true>;
