//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MooseVariableFV.h"
#include <typeinfo>
#include "TimeIntegrator.h"
#include "NonlinearSystemBase.h"
#include "DisplacedSystem.h"
#include "Assembly.h"
#include "MooseVariableData.h"

template <typename OutputType>
MooseVariableFV<OutputType>::MooseVariableFV(const InputParameters & parameters)
  : MooseVariableFVBase(parameters)
{
  _element_data = libmesh_make_unique<MooseVariableData<OutputType>>(*this,
                                                                     _sys,
                                                                     _tid,
                                                                     Moose::ElementType::Element,
                                                                     _assembly.qRule(),
                                                                     _assembly.qRuleFace(),
                                                                     _assembly.node(),
                                                                     _assembly.elem());
  _neighbor_data =
      libmesh_make_unique<MooseVariableData<OutputType>>(*this,
                                                         _sys,
                                                         _tid,
                                                         Moose::ElementType::Neighbor,
                                                         _assembly.qRuleNeighbor(), // Place holder
                                                         _assembly.qRuleNeighbor(),
                                                         _assembly.nodeNeighbor(),
                                                         _assembly.neighbor());
}

template <typename OutputType>
const std::set<SubdomainID> &
MooseVariableFV<OutputType>::activeSubdomains() const
{
  return _sys.system().variable(_var_num).active_subdomains();
}

template <typename OutputType>
Moose::VarFieldType
MooseVariableFV<OutputType>::fieldType() const
{
  if (std::is_same<OutputType, Real>::value)
    return Moose::VarFieldType::VAR_FIELD_STANDARD;
  else if (std::is_same<OutputType, RealVectorValue>::value)
    return Moose::VarFieldType::VAR_FIELD_VECTOR;
  else if (std::is_same<OutputType, RealEigenVector>::value)
    return Moose::VarFieldType::VAR_FIELD_ARRAY;
  else
    mooseError("Unknown variable field type");
}

template <typename OutputType>
bool
MooseVariableFV<OutputType>::activeOnSubdomain(SubdomainID subdomain) const
{
  return _sys.system().variable(_var_num).active_on_subdomain(subdomain);
}

template <typename OutputType>
void
MooseVariableFV<OutputType>::clearDofIndices()
{
  _element_data->clearDofIndices();
}

template <typename OutputType>
void
MooseVariableFV<OutputType>::prepare()
{
  _element_data->prepare();
}
template <typename OutputType>
void
MooseVariableFV<OutputType>::prepareFace()
{
  _element_data->prepare();
  _neighbor_data->prepare();
}

template <typename OutputType>
void
MooseVariableFV<OutputType>::getDofIndices(const Elem * elem,
                                           std::vector<dof_id_type> & dof_indices) const
{
  _element_data->getDofIndices(elem, dof_indices);
}

template <typename OutputType>
typename MooseVariableFV<OutputType>::OutputData
MooseVariableFV<OutputType>::getElementalValue(const Elem * elem, unsigned int idx) const
{
  return _element_data->getElementalValue(elem, Moose::Current, idx);
}

template <typename OutputType>
typename MooseVariableFV<OutputType>::OutputData
MooseVariableFV<OutputType>::getElementalValueOld(const Elem * elem, unsigned int idx) const
{
  return _element_data->getElementalValue(elem, Moose::Old, idx);
}

template <typename OutputType>
typename MooseVariableFV<OutputType>::OutputData
MooseVariableFV<OutputType>::getElementalValueOlder(const Elem * elem, unsigned int idx) const
{
  return _element_data->getElementalValue(elem, Moose::Older, idx);
}

template <typename OutputType>
void
MooseVariableFV<OutputType>::insert(NumericVector<Number> & residual)
{
  _element_data->insert(residual);
}

template <typename OutputType>
void
MooseVariableFV<OutputType>::add(NumericVector<Number> & residual)
{
  _element_data->add(residual);
}

template <typename OutputType>
void
MooseVariableFV<OutputType>::addSolution(const DenseVector<Number> & v)
{
  _element_data->addSolution(_sys.solution(), v);
}

template <typename OutputType>
void
MooseVariableFV<OutputType>::addSolutionNeighbor(const DenseVector<Number> & v)
{
  _neighbor_data->addSolution(_sys.solution(), v);
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::DoFValue &
MooseVariableFV<OutputType>::dofValues()
{
  return _element_data->dofValues();
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::DoFValue &
MooseVariableFV<OutputType>::dofValuesOld()
{
  return _element_data->dofValuesOld();
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::DoFValue &
MooseVariableFV<OutputType>::dofValuesOlder()
{
  return _element_data->dofValuesOlder();
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::DoFValue &
MooseVariableFV<OutputType>::dofValuesPreviousNL()
{
  return _element_data->dofValuesPreviousNL();
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::DoFValue &
MooseVariableFV<OutputType>::dofValuesNeighbor()
{
  return _neighbor_data->dofValues();
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::DoFValue &
MooseVariableFV<OutputType>::dofValuesOldNeighbor()
{
  return _neighbor_data->dofValuesOld();
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::DoFValue &
MooseVariableFV<OutputType>::dofValuesOlderNeighbor()
{
  return _neighbor_data->dofValuesOlder();
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::DoFValue &
MooseVariableFV<OutputType>::dofValuesPreviousNLNeighbor()
{
  return _neighbor_data->dofValuesPreviousNL();
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::DoFValue &
MooseVariableFV<OutputType>::dofValuesDot()
{
  return _element_data->dofValuesDot();
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::DoFValue &
MooseVariableFV<OutputType>::dofValuesDotDot()
{
  return _element_data->dofValuesDotDot();
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::DoFValue &
MooseVariableFV<OutputType>::dofValuesDotOld()
{
  return _element_data->dofValuesDotOld();
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::DoFValue &
MooseVariableFV<OutputType>::dofValuesDotDotOld()
{
  return _element_data->dofValuesDotDotOld();
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::DoFValue &
MooseVariableFV<OutputType>::dofValuesDotNeighbor()
{
  return _neighbor_data->dofValuesDot();
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::DoFValue &
MooseVariableFV<OutputType>::dofValuesDotDotNeighbor()
{
  return _neighbor_data->dofValuesDotDot();
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::DoFValue &
MooseVariableFV<OutputType>::dofValuesDotOldNeighbor()
{
  return _neighbor_data->dofValuesDotOld();
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::DoFValue &
MooseVariableFV<OutputType>::dofValuesDotDotOldNeighbor()
{
  return _neighbor_data->dofValuesDotDotOld();
}

template <typename OutputType>
const MooseArray<Number> &
MooseVariableFV<OutputType>::dofValuesDuDotDu()
{
  return _element_data->dofValuesDuDotDu();
}

template <typename OutputType>
const MooseArray<Number> &
MooseVariableFV<OutputType>::dofValuesDuDotDotDu()
{
  return _element_data->dofValuesDuDotDotDu();
}

template <typename OutputType>
const MooseArray<Number> &
MooseVariableFV<OutputType>::dofValuesDuDotDuNeighbor()
{
  return _neighbor_data->dofValuesDuDotDu();
}

template <typename OutputType>
const MooseArray<Number> &
MooseVariableFV<OutputType>::dofValuesDuDotDotDuNeighbor()
{
  return _neighbor_data->dofValuesDuDotDotDu();
}

template <typename OutputType>
void
MooseVariableFV<OutputType>::prepareIC()
{
  _element_data->prepareIC();
}

template <typename OutputType>
void
MooseVariableFV<OutputType>::computeElemValues()
{
  _element_data->setGeometry(Moose::Volume);
  _element_data->computeValues();
}

template <typename OutputType>
void
MooseVariableFV<OutputType>::computeElemValuesFace()
{
  _element_data->setGeometry(Moose::Face);
  _element_data->computeValues();
}

template <typename OutputType>
void
MooseVariableFV<OutputType>::computeNeighborValuesFace()
{
  _neighbor_data->setGeometry(Moose::Face);
  _neighbor_data->computeValues();
}

template <typename OutputType>
void
MooseVariableFV<OutputType>::computeNeighborValues()
{
  _neighbor_data->setGeometry(Moose::Volume);
  _neighbor_data->computeValues();
}

template <typename OutputType>
OutputType
MooseVariableFV<OutputType>::getValue(const Elem * elem) const
{
  std::vector<dof_id_type> dof_indices;
  _dof_map.dof_indices(elem, dof_indices, _var_num);
  mooseAssert(dof_indices.size() == 1, "Wrong size for dof indices");
  OutputType value = (*_sys.currentSolution())(dof_indices[0]);
  return value;
}

template <>
RealEigenVector
MooseVariableFV<RealEigenVector>::getValue(const Elem * elem) const
{
  std::vector<dof_id_type> dof_indices;
  _dof_map.dof_indices(elem, dof_indices, _var_num);

  RealEigenVector value(_count);
  mooseAssert(dof_indices.size() == 1, "Wrong size for dof indices");
  unsigned int n = 0;
  for (unsigned int j = 0; j < _count; j++)
  {
    value(j) = (*_sys.currentSolution())(dof_indices[0] + n);
    n += _dof_indices.size();
  }

  return value;
}

template <typename OutputType>
typename OutputTools<OutputType>::OutputGradient
MooseVariableFV<OutputType>::getGradient(const Elem * elem) const
{
  return 0;
}

template <>
RealVectorArrayValue
MooseVariableFV<RealEigenVector>::getGradient(
    const Elem * elem, const std::vector<std::vector<RealVectorValue>> & grad_phi) const
{
  std::vector<dof_id_type> dof_indices;
  _dof_map.dof_indices(elem, dof_indices, _var_num);

  RealVectorArrayValue value(_count, LIBMESH_DIM);
  if (isNodal())
  {
    for (unsigned int i = 0; i < dof_indices.size(); ++i)
      for (unsigned int j = 0; j < _count; ++j)
        for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
        {
          // The zero index is because we only have one point that the phis are evaluated at
          value(j, k) += grad_phi[i][0](k) * (*_sys.currentSolution())(dof_indices[i] + j);
        }
  }
  else
  {
    mooseAssert(dof_indices.size() == 1, "Wrong size for dof indices");
  }

  return value;
}

template <typename OutputType>
void
MooseVariableFV<OutputType>::setDofValue(const OutputData & value, unsigned int index)
{
  _element_data->setDofValue(value, index);
}

template <typename OutputType>
void
MooseVariableFV<OutputType>::setDofValues(const DenseVector<OutputData> & values)
{
  _element_data->setDofValues(values);
}

template <typename OutputType>
bool
MooseVariableFV<OutputType>::isVector() const
{
  return std::is_same<OutputType, RealVectorValue>::value;
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::FieldVariablePhiSecond &
MooseVariableFV<OutputType>::secondPhi() const
{
  return _element_data->secondPhi();
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::FieldVariablePhiCurl &
MooseVariableFV<OutputType>::curlPhi() const
{
  return _element_data->curlPhi();
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::FieldVariablePhiSecond &
MooseVariableFV<OutputType>::secondPhiFace() const
{
  return _element_data->secondPhiFace();
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::FieldVariablePhiCurl &
MooseVariableFV<OutputType>::curlPhiFace() const
{
  return _element_data->curlPhiFace();
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::FieldVariablePhiSecond &
MooseVariableFV<OutputType>::secondPhiNeighbor() const
{
  return _neighbor_data->secondPhi();
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::FieldVariablePhiCurl &
MooseVariableFV<OutputType>::curlPhiNeighbor() const
{
  return _neighbor_data->curlPhi();
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::FieldVariablePhiSecond &
MooseVariableFV<OutputType>::secondPhiFaceNeighbor() const
{
  return _neighbor_data->secondPhiFace();
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::FieldVariablePhiCurl &
MooseVariableFV<OutputType>::curlPhiFaceNeighbor() const
{
  return _neighbor_data->curlPhiFace();
}

template <typename OutputType>
bool
MooseVariableFV<OutputType>::usesSecondPhi() const
{
  return _element_data->usesSecondPhi();
}

template <typename OutputType>
bool
MooseVariableFV<OutputType>::usesSecondPhiNeighbor() const
{
  return _neighbor_data->usesSecondPhi();
}

template <typename OutputType>
bool
MooseVariableFV<OutputType>::computingCurl() const
{
  return _element_data->computingCurl();
}

template class MooseVariableFV<Real>;
template class MooseVariableFV<RealVectorValue>;
template class MooseVariableFV<RealEigenVector>;
