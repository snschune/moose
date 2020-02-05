#pragma once

#include "FVFluxKernel.h"

class FVDiffusion : public FVFluxKernel
{
public:
  FVDiffusion(const InputParameters & params) : FVFluxKernel(params),  _coeff_left(getMaterialProperty("coeff")), _coeff_right(getNeighborMaterialProperty("coeff"){};

protected:
  virtual Real computeQpResidual() override
  {
    Real k = (_coeff_left[_qp] + _coeff_right[_qp]) / 2;
    return _normal[_qp] * k * _grad_u_face[_qp];
  }

  const Real & _coeff_left;
  const Real & _coeff_right;
};
