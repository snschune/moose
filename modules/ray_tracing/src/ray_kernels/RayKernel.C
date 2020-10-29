//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RayKernel.h"

// Local includes
#include "RayTracingStudy.h"

// MOOSE includes
#include "Assembly.h"
#include "NonlinearSystemBase.h"

template <typename T>
InputParameters
RayKernelTempl<T>::validParams()
{
  auto params = IntegralRayKernelBase::validParams();
  params += TaggingInterface::validParams();

  params.template addRequiredParam<NonlinearVariableName>(
      "variable", "The name of the variable that this RayKernel operates on");

  return params;
}

template <typename T>
RayKernelTempl<T>::RayKernelTempl(const InputParameters & params)
  : IntegralRayKernelBase(params),
    MooseVariableInterface<T>(this,
                              false,
                              "variable",
                              Moose::VarKindType::VAR_NONLINEAR,
                              std::is_same<T, Real>::value ? Moose::VarFieldType::VAR_FIELD_STANDARD
                                                           : Moose::VarFieldType::VAR_FIELD_VECTOR),
    TaggingInterface(this),
    _var(*this->mooseVariable()),
    _u(_is_implicit ? _var.sln() : _var.slnOld()),
    _grad_u(_is_implicit ? _var.gradSln() : _var.gradSlnOld()),
    _test(_var.phi()),
    _grad_test(_var.gradPhi()),
    _phi(_assembly.phi(_var)),
    _grad_phi(_assembly.gradPhi(_var))
{
  // We do not allow RZ/RSPHEICAL because in the context of these coord
  // systems there is no way to represent a line source - we would end up
  // with a plane/surface source or a volumetric source, respectively.
  // This is also why we do not multiply by _coord[_qp] in any of the
  // integrations that follow.
  for (const auto & subdomain_id : _mesh.meshSubdomains())
    if (_fe_problem.getCoordSystem(subdomain_id) != Moose::COORD_XYZ)
      mooseError(_error_prefix, ": Not valid on coordinate systems other than XYZ");

  addMooseVariableDependency(this->mooseVariable());
}

template <typename T>
void
RayKernelTempl<T>::computeIntegral()
{
  mooseAssert(_current_subdomain_id == _assembly.currentSubdomainID(), "Subdomain IDs not in sync");

  if (_study.computingJacobian())
    computeJacobian();
  else if (_study.computingResidual())
    computeResidual();
  else
    mooseError(_error_prefix, " Should not be computing outside of residual/Jacobian");
}

template <typename T>
void
RayKernelTempl<T>::computeResidual()
{
  prepareVectorTag(_assembly, _var.number());

  precalculateResidual();
  for (_i = 0; _i < _test.size(); ++_i)
    for (_qp = 0; _qp < _q_point.size(); ++_qp)
      _local_re(_i) += _JxW[_qp] * computeQpResidual();

  accumulateTaggedLocalResidual();
}

template <typename T>
void
RayKernelTempl<T>::computeJacobian()
{
  if (!isImplicit())
    return;

  _subproblem.prepareShapes(_var.number(), _tid);

  precalculateJacobian();

  const auto & ce = _fe_problem.couplingEntries(_tid);
  for (const auto & it : ce)
  {
    MooseVariableFEBase & ivariable = *(it.first);
    MooseVariableFEBase & jvariable = *(it.second);

    const auto ivar = ivariable.number();
    if (ivar != _var.number())
      continue;

    if (!jvariable.activeOnSubdomain(_current_subdomain_id))
      continue;

    if (jvariable.count() > 1)
      mooseError(_error_prefix, ": ArrayMooseVariable objects are not coupleable");

    const auto jvar = jvariable.number();

    prepareMatrixTag(_assembly, _var.number(), jvar);
    if (_local_ke.m() != _test.size())
      return;

    if (ivar == jvar)
      precalculateJacobian();
    else
      precalculateOffDiagJacobian(jvar);

    for (_i = 0; _i < _test.size(); _i++)
      for (_j = 0; _j < _phi.size(); _j++)
        for (_qp = 0; _qp < _JxW.size(); _qp++)
          _local_ke(_i, _j) += _JxW[_qp] * computeQpOffDiagJacobian(jvar);

    accumulateTaggedLocalMatrix();
  }
}

template class RayKernelTempl<Real>;
// Not implementing this until there is a use case and tests for it!
// template class RayKernelTempl<RealVectorValue>;
