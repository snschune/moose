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
#ifndef PicardFunctionDT_H
#define PicardFunctionDT_H

#include "FunctionDT.h"

class PicardFunctionDT;

template<>
InputParameters validParams<PicardFunctionDT>();

/**
 * FunctionDT that also detects Picard Divergence
 */
class PicardFunctionDT : public FunctionDT
{
public:
  PicardFunctionDT(const InputParameters & parameters);

  /**
   * If the time step converged
   * @return true if converged, otherwise false
   */
  virtual bool converged() { return _converged && !picardDivergence(); }

protected:
  /// detects the failure of the Picard iteration
  bool picardDivergence();

  /// the maximum growth factor permitted before time step cutback w.r.t. the inital picard residual
  Real _max_residual_growth_factor;
};

#endif /* PicardFunctionDT_H */
