//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADRayKernel.h"

// Local includes
#include "RayTracingStudy.h"

// MOOSE includes
#include "Assembly.h"
#include "NonlinearSystemBase.h"
#include "ADUtils.h"

template <typename T>
InputParameters
ADRayKernelTempl<T>::validParams()
{
  auto params = IntegralRayKernelBase::validParams();
  params += TaggingInterface::validParams();

  params.addRequiredParam<NonlinearVariableName>(
      "variable", "The name of the variable that this ADRayKernel operates on");

  return params;
}

template <typename T>
ADRayKernelTempl<T>::ADRayKernelTempl(const InputParameters & params)
  : IntegralRayKernelBase(params),
    MooseVariableInterface<T>(this,
                              false,
                              "variable",
                              Moose::VarKindType::VAR_NONLINEAR,
                              std::is_same<T, Real>::value ? Moose::VarFieldType::VAR_FIELD_STANDARD
                                                           : Moose::VarFieldType::VAR_FIELD_VECTOR),
    TaggingInterface(this),
    _var(*this->mooseVariable()),
    _test(_var.phi()),
    _grad_test(_var.adGradPhi()),
    _u(_var.adSln()),
    _grad_u(_var.adGradSln()),
    _phi(_assembly.phi(_var)),
    _grad_phi(_assembly.template adGradPhi<T>(_var))
{
  // We do not allow RZ/RSPHEICAL because in the context of these coord
  // systems there is no way to represent a line source - we would end up
  // with a plane/surface source or a volumetric source, respectively.
  // This is also why we do not multiply by _coord[_qp] in any of the
  // integrations that follow.
  for (const auto & subdomain_id : _mesh.meshSubdomains())
    if (_fe_problem.getCoordSystem(subdomain_id) != Moose::COORD_XYZ)
      mooseError(_error_prefix, ": Not valid on coordinate systems other than XYZ");

  _subproblem.haveADObjects(true);

  addMooseVariableDependency(this->mooseVariable());
}

template <typename T>
void
ADRayKernelTempl<T>::computeIntegral()
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
ADRayKernelTempl<T>::computeResidual()
{
  prepareVectorTag(_assembly, _var.number());

  precalculateResidual();
  for (_qp = 0; _qp < _JxW.size(); _qp++)
    for (_i = 0; _i < _test.size(); _i++)
      _local_re(_i) += raw_value(_JxW[_qp] * computeQpResidual());

  accumulateTaggedLocalResidual();
}

template <typename T>
void
ADRayKernelTempl<T>::computeJacobian()
{
  if (!isImplicit())
    return;

  _subproblem.prepareShapes(_var.number(), _tid);

  _residuals.resize(_test.size(), 0);
  for (auto & r : _residuals)
    r = 0;

  precalculateResidual();
  for (_qp = 0; _qp < _JxW.size(); _qp++)
    for (_i = 0; _i < _test.size(); _i++)
      _residuals[_i] += _JxW[_qp] * computeQpResidual();

  const auto & ce = _assembly.couplingEntries();
  for (const auto & it : ce)
  {
    MooseVariableFEBase & ivariable = *(it.first);
    MooseVariableFEBase & jvariable = *(it.second);

    const auto ivar = ivariable.number();
    if (ivar != _var.number())
      continue;

    if (!jvariable.activeOnSubdomain(_current_subdomain_id))
      continue;

    const auto jvar = jvariable.number();
    const auto ad_offset =
        Moose::adOffset(jvar, _nl.getMaxVarNDofsPerElem(), Moose::ElementType::Element);

    prepareMatrixTag(_assembly, ivar, jvar);

    if (_local_ke.m() != _test.size() || _local_ke.n() != jvariable.phiSize())
      continue;

    for (_i = 0; _i < _test.size(); _i++)
      for (_j = 0; _j < jvariable.phiSize(); _j++)
        _local_ke(_i, _j) += _residuals[_i].derivatives()[ad_offset + _j];

    accumulateTaggedLocalMatrix();
  }
}

template class ADRayKernelTempl<Real>;

// Not implementing this until there is a use case and tests for it!
// template class ADRayKernelTempl<RealVectorValue>;
