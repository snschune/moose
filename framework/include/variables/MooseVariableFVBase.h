//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseVariableBase.h"

namespace libMesh
{
template <typename>
class NumericVector;
}

class FaceInfo;

class MooseVariableFVBase;
template <>
InputParameters validParams<MooseVariableFVBase>();

class MooseVariableFVBase : public MooseVariableBase
{
public:
  static InputParameters validParams();

  MooseVariableFVBase(const InputParameters & parameters);

  /// Clear out the dof indices.  We do this in case this variable is not going to be prepared.
  virtual void clearDofIndices() = 0;

  /// Prepare the elemental degrees of freedom
  virtual void prepare() = 0;

  /// Prepare the element+neighbor element degrees of freedom for the given face
  virtual void prepareFace(const FaceInfo & fi) = 0;

  /// Prepare the initial condition
  virtual void prepareIC() = 0;

  /// Filed type of this variable
  virtual Moose::VarFieldType fieldType() const = 0;

  /// returns true if this is a vector-valued element, false otherwise.
  virtual bool isVector() const = 0;

  virtual bool isNodal() const override { return false; }

  /// Current element this variable is being evaluated at in a volumetric/elemental
  /// or face/flux loop.
  virtual const Elem * const & currentElem() const = 0;
  /// Current neighbor element this variable is being evaluated at in a
  /// face/flux loop - this is nullptr for an elemental loop.
  virtual const Elem * const & currentNeighbor() const = 0;

  /// The subdomains the variable is active on
  virtual const std::set<SubdomainID> & activeSubdomains() const = 0;
  /// returns true if the variable is active on the given subdomain
  virtual bool activeOnSubdomain(SubdomainID subdomain) const = 0;

  /// Compute values at interior quadrature points
  virtual void computeElemValues() = 0;
  /// Compute values at face quadrature points for the element+neighbor (both
  /// sides of the face).
  virtual void computeFaceValues(const FaceInfo & fi) = 0;

  virtual const std::vector<dof_id_type> & dofIndicesNeighbor() const = 0;

  virtual void insert(NumericVector<Number> & residual) = 0;
  virtual void add(NumericVector<Number> & residual) = 0;
};
