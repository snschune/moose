/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

// MOOSE includes
#include "PicardFunctionDT.h"
#include "Transient.h"

template<>
InputParameters validParams<PicardFunctionDT>()
{
  InputParameters params = validParams<FunctionDT>();
  params.addClassDescription("FunctionDT that also detects Picard divergence");
  params.addRangeCheckedParam<Real>("max_residual_growth_factor", std::numeric_limits<Real>::max(), "max_residual_growth_factor>1", "The maximum ratio between any Picard residual and the initial Picard residual before timestep is cut."
                                    " This catches diverging Picard solutions");
  return params;
}

PicardFunctionDT::PicardFunctionDT(const InputParameters & parameters) :
    FunctionDT(parameters),
    _max_residual_growth_factor(getParam<Real>("max_residual_growth_factor"))
{
}

bool
PicardFunctionDT::picardDivergence()
{
  if (_executioner.numPicardIts() == 1)
    return false;

  return _executioner.picardResidualDrop() > _max_residual_growth_factor;
}
