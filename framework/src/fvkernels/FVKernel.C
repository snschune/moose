
#include "FVKernel.h"
#include "Assembly.h"

#include "ComputeFVFaceResidualsThread.h"

FVFluxKernel::FVFluxKernel(const InputParameters & params)
  : MooseObject(params),
    TaggingInterface(this),
    TransientInterface(this),
    BlockRestrictable(this),
    TwoMaterialPropertyInterface(this, blockIDs(), {}),
    NeighborCoupleableMooseVariableDependencyIntermediateInterface(this, false, false),
    NeighborMooseVariableInterface(
        this, false, Moose::VarKindType::VAR_NONLINEAR, Moose::VarFieldType::VAR_FIELD_STANDARD),
    _var(*mooseVariable()),
    _tid(params.get<THREAD_ID>("_tid")),
    _assembly(_subproblem.assembly(_tid)),
    _u_left(_is_implicit ? _var.sln() : _var.slnOld()),
    _u_right(_is_implicit ? _var.slnNeighbor() : _var.slnOldNeighbor()),
    _grad_u_left(_is_implicit ? _var.gradSln() : _var.gradSlnOld()),
    _grad_u_right(_is_implicit ? _var.gradSlnNeighbor() : _var.gradSlnOldNeighbor()),
    _normal(_assembly.normals())
{
}

Real
FVFluxKernel::computeResidual(const FaceInfo & fi)
{
  _face_info = &fi;
  auto r = fi.faceArea() * computeQpResidual(fi);

  if (ownLeftElem())
  {
    prepareVectorTag(_assembly, _var.number());
    _local_re(0) = r;
    accumulateTaggedLocalResidual();
  }
  if (ownRightElem())
  {
    prepareVectorTagNeighbor(_assembly, _var.number());
    _local_re(0) = -r;
    accumulateTaggedLocalResidual();
  }
}
