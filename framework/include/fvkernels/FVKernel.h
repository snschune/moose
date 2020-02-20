#pragma once

#include "MooseObject.h"
#include "TaggingInterface.h"
#include "TransientInterface.h"
#include "BlockRestrictable.h"
#include "NeighborCoupleable.h"
#include "TwoMaterialPropertyInterface.h"
#include "NeighborMooseVariableInterface.h"
#include "MooseVariableDependencyInterface.h"

class FaceInfo;

// TODO: implement this eventually
class FVBoundaryCondition
{
};

class FVKernel : public MooseObject,
                 public TaggingInterface,
                 public TransientInterface,
                 public BlockRestrictable,
                 public MooseVariableDependencyInterface
{
public:
  FVKernel(const InputParameters & params);
};

class FVFluxKernel : public FVKernel,
                     public TwoMaterialPropertyInterface,
                     public NeighborCoupleable,
                     public NeighborMooseVariableInterface<Real>
{
public:
  FVFluxKernel(const InputParameters & params);

  virtual void computeResidual(const FaceInfo & fi);

protected:
  // material properties will be initialized on the face.  Reconstructed
  // solutions will have been performed previous to this call and all coupled
  // variables and _u, _grad_u, etc. will have their reconstructed values
  // extrapolated to/at the face.  Material properties will also have been
  // computed on the face using the face-reconstructed variable values.
  //
  virtual Real computeQpResidual(const FaceInfo & fi) = 0;

  MooseVariable & _var;
  THREAD_ID _tid;
  Assembly & _assembly;

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
