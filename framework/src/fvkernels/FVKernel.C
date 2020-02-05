
#include "FVKernel.h"

FVKernel::FVKernel(const InputParameters & params)
  : MooseObject(params), TaggingInterface(this), TransientInterface(this)
{
}

FVFluxKernel::FVFluxKernel(const InputParameters & params)
  : MooseObject(params),
    TaggingInterface(this),
    TransientInterface(this),
    _var(*mooseVariable()),
    _u_left(_is_implicit ? _var.sln() : _var.slnOld()),
    _u_right(_is_implicit ? _var.slnNeighbor() : _var.slnOldNeighbor()),
    _grad_u_left(_is_implicit ? _var.gradSln() : _var.gradSlnOld()),
    _grad_u_right(_is_implicit ? _var.gradSlnNeighbor() : _var.gradSlnOldNeighbor()),
    _normal(_assembly.normals()[0]){};

Real
FVFluxKernel::computeResidual(const FaceInfo & fi)
{
  _face_info = &fi;
  auto r = _face_area * computeQpResidual();

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
