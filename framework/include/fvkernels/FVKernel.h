#pragma once

class FVKernel : public MooseObject, public TaggingInterface
{
public:
  FVKernel(const InputParameters & params);
  ADReal computeResidual() = 0;
};

class ConvectiveFluxAverage : public FVInterpMethod
{
  Real interpolate(FVKernelFace & k)
  {
    k.initLeft();
    auto v_left = k.oneSidedConvectiveFlux();
    k.initRight();
    auto v_right = k.oneSidedConvectiveFlux();
    return (v_left + v_right) / 2;
  }
};

class UpwindMethod : public FVInterpMethod
{
  Real interpolate(FVKernelFace & k)
  {
    auto cv = k.convectiveFlux();
    if (cv > 0)
      k.initLeft();
    else
      k.initRight;
    return k.computeQpResidual();
  }
};

class LinearUpwind : public FVInterpMethod
{
  Real interpolate(FVKernelFace k)
  {
    k.initLeft();
    auto r1 = k.computeQpResidual();
    k.initRight();
    auto r2 = k.computeQpResidual();

    Real d1 = k.left();
    Real d2 =
        // ... get convective flux, etc.

        return d1 / (d1 + d2) * r2 + r1 // ...
  }
};

class ResidualAverage : public FVInterpMethod
{
  Real interpolate(FVKernelFace & k)
  {
    k.initLeft();
    auto v_left = k.computeQpResidual();
    k.initRight();
    auto v_right = k.computeQpResidual();
    return (v_left + v_right) / 2;
  }
};

// TODO: need to implement initLeft and initRight functionality - this should
// control which material property and coupled variable values are exposed via
// _u, _grad_u, _random_mat_prop_name, etc. - this will require writing some
// sort of new FVCoupleable and FVMaterialProperty interfaces that toggle
// between left and right elements' info.
class FVFluxKernel : public FVKernel
{
public:
  FVFluxKernel(const InputParameters & params)
    : _residual_interp_method(getParam(...)), _matprop_iface(this, {}, {}){};

  Real computeResidual()
  {
    auto r = _face_area * computeQpResidual();

    if (ownElement())
    {
      prepareVectorTag(_assembly, _var.number());
      _local_re(0) = r;
      accumulateTaggedLocalResidual();
    }
    if (ownNeighbor())
    {
      prepareVectorTagNeighbor(_assembly, _var.number());
      _local_re(0) = -r;
      accumulateTaggedLocalResidual();
    }
  }

  Real convectiveFlux() { return _convective_flux; }

protected:
  virtual Real oneSidedConvectiveFlux()
  {
    return 0;
    // i.e. return _normal * _rho * _vel;
    // or return _normal * _vel;
    // or return _normal * _fluid_prop.rho_from_p_T(_pressure, _temperature) * _vel;
    // etc.
  }

  // material properties will be initialized on the face.  Reconstructed
  // solutions will have been performed previous to this call and all coupled
  // variables and _u, _grad_u, etc. will have their reconstructed values
  // extrapolated to/at the face. This is possible because of the
  // oneSidedConvectiveFlux function that helps us interpolate the advective
  // flux quantity at the face so we can use it to reconstruct all these
  // solutions.
  virtual Real computeQpResidual()
  {
    // for convective terms, e.g.: return convectiveFlux() * _u;
    // for other terms, e.g.: -1 * _normal *  _matprop * _grad_u
  }

private:
  Real _convective_flux = 0;

  FVInterpMethod & _residual_interp_method;
  FVInterpMethod * _flux_interp_method = nullptr;
};
