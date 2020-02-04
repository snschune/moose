#pragma once

class FVKernel : public MooseObject, public TaggingInterface, public TransientInterface
{
public:
  FVKernel(const InputParameters & params);
  virtual Real computeResidual(const FaceInfo & fi) = 0;

protected:
  virtual Real computeQpResidual() = 0;
};

class FVFluxKernel : public FVKernel, public NeighborCoupleable, public TwoMaterialPropertyInterface
{
public:
  FVFluxKernel(const InputParameters & params);

  virtual Real computeResidual(const FaceInfo & fi);

protected:
  // material properties will be initialized on the face.  Reconstructed
  // solutions will have been performed previous to this call and all coupled
  // variables and _u, _grad_u, etc. will have their reconstructed values
  // extrapolated to/at the face.  Material properties will also have been
  // computed on the face using the face-reconstructed variable values.
  //
  //   virtual Real computeQpResidual() override { ///... }

  const FaceInfo * _face_info;
  MooseVariable & _var;
  const VariableValue & _u_left;
  const VariableValue & _u_right;
  const VariableGradient & _grad_u_left;
  const VariableGradient & _grad_u_right;
  const VariableGradient _grad_u_face;

private:
  bool ownElement()
  {
    // returns true if this processor owns the (left) element.
    // TODO: implement this
    return true;
  }
  bool ownNeighbor()
  {
    // returns true if this processor owns the (right) neighbor.
    // TODO: implement this
    return true;
  }
};
