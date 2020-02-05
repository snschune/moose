#pragma once

class FVFluxKernel : public MooseObject,
                     public TaggingInterface,
                     public TransientInterface,
                     public NeighborCoupleable,
                     public TwoMaterialPropertyInterface
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
  virtual Real computeQpResidual(const FaceInfo & fi) = 0;

  MooseVariable & _var;

  const VariableValue & _u_left;
  const VariableValue & _u_right;
  const VariableGradient & _grad_u_left;
  const VariableGradient & _grad_u_right;
  const MooseArray<Point> & _normal;

  const FaceInfo * _face_info = nullptr;

private:
  bool ownLeftElem()
  {
    // returns true if this processor owns the (left) element.
    // TODO: implement this
    return true;
  }
  bool ownRightElem()
  {
    // returns true if this processor owns the (right) neighbor.
    // TODO: implement this
    return true;
  }
};
