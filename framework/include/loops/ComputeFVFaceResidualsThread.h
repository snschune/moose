//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ParallelUniqueId.h"
#include "MooseMesh.h"
#include "MooseTypes.h"
#include "MooseException.h"
#include "libmesh/libmesh_exceptions.h"
#include "libmesh/elem.h"

class FaceInfo
{
public:
  FaceInfo() {}
  Real faceArea();
  const RealVectorValue & unitNormalVec();
  const std::vector<Point> & faceNodes();

  const Point & faceCentroid();
  const Point & leftCentroid();
  const Point & rightCentroid();

  const Elem & leftElem();
  const Elem & rightElem();
  unsigned int leftSide();
  unsigned int rightSide();
};

/**
 * Base class for assembly-like calculations.
 */
template <typename RangeType>
class ComputeFVFaceResidualsThread
{
public:
  ComputeFVFaceResidualsThread(MooseMesh & mesh);

  ComputeFVFaceResidualsThread(ComputeFVFaceResidualsThread & x, Threads::split split);

  virtual ~ComputeFVFaceResidualsThread();

  virtual void operator()(const RangeType & range, bool bypass_threading = false);

  /**
   * Called before the element range loop
   */
  virtual void pre();

  /**
   * Called after the element range loop
   */
  virtual void post();

  /**
   * Called before the boundary assembly
   *
   * @param elem - The element we are checking is on the boundary.
   * @param side - The side of the element in question.
   * @param bnd_id - ID of the boundary we are at
   */
  virtual void onFace(const FaceInfo & fi);

  /**
   * Called before the boundary assembly
   *
   * @param bnd_id - ID of the boundary we are at
   */
  virtual void onBoundary(const FaceInfo & fi, BoundaryID boundary);

  /**
   * Called every time the current subdomain changes (i.e. the subdomain of _this_ element
   * is not the same as the subdomain of the last element).  Beware of over-using this!
   * You might think that you can do some expensive stuff in here and get away with it...
   * but there are applications that have TONS of subdomains....
   */
  virtual void subdomainChanged();

  /**
   * Called every time the neighbor subdomain changes (i.e. the subdomain of _this_ neighbor
   * is not the same as the subdomain of the last neighbor).  Beware of over-using this!
   * You might think that you can do some expensive stuff in here and get away with it...
   * but there are applications that have TONS of subdomains....
   */
  virtual void neighborSubdomainChanged(){};

  /**
   * Called if a MooseException is caught anywhere during the computation.
   * The single input parameter taken is a MooseException object.
   */
  virtual void caughtMooseException(MooseException &){};

protected:
  MooseMesh & _mesh;
  THREAD_ID _tid;

  /// The subdomain for the current element
  SubdomainID _subdomain;

  /// The subdomain for the last element
  SubdomainID _old_subdomain;

  /// The subdomain for the current neighbor
  SubdomainID _neighbor_subdomain;

  /// The subdomain for the last neighbor
  SubdomainID _old_neighbor_subdomain;
};

template <typename RangeType>
ComputeFVFaceResidualsThread<RangeType>::ComputeFVFaceResidualsThread(MooseMesh & mesh)
  : _mesh(mesh)
{
}

template <typename RangeType>
ComputeFVFaceResidualsThread<RangeType>::ComputeFVFaceResidualsThread(
    ComputeFVFaceResidualsThread & x, Threads::split /*split*/)
  : _mesh(x._mesh)
{
}

template <typename RangeType>
ComputeFVFaceResidualsThread<RangeType>::~ComputeFVFaceResidualsThread()
{
}

// the vector<faceinfo> data structure needs to be built such that for all
// sides on an interface between two subdomains, the elements of the same
// subdomain are used consistently for all the "elem" (i.e. not "neighbor")
// parameters in order to avoid jumping back and forth along the boundary
// between using one or the other subdomains' FV kernels unpredictably.
template <typename RangeType>
void
ComputeFVFaceResidualsThread<RangeType>::operator()(const RangeType & range, bool bypass_threading)
{
  try
  {
    try
    {
      ParallelUniqueId puid;
      _tid = bypass_threading ? 0 : puid.id;

      pre();

      _subdomain = Moose::INVALID_BLOCK_ID;
      _neighbor_subdomain = Moose::INVALID_BLOCK_ID;

      typename RangeType::const_iterator faceinfo = range.begin();
      for (faceinfo = range.begin(); faceinfo != range.end(); ++faceinfo)
      {
        const Elem & elem = faceinfo->leftElem();
        unsigned int side = faceinfo->leftSide();

        _old_subdomain = _subdomain;
        _subdomain = elem->subdomain_id();
        if (_subdomain != _old_subdomain)
          subdomainChanged();

        // get elem's face residual contribution to it's neighbor
        onFace(*faceinfo);

        // boundary faces only border one element and so only contribute to
        // one element's residual
        std::vector<BoundaryID> boundary_ids = _mesh.getBoundaryIDs(&elem, side);
        for (auto it = boundary_ids.begin(); it != boundary_ids.end(); ++it)
          onBoundary(*faceinfo, *it);
      } // range

      post();
    }
    catch (libMesh::LogicError & e)
    {
      throw MooseException("We caught a libMesh error");
    }
  }
  catch (MooseException & e)
  {
    caughtMooseException(e);
  }
}

template <typename RangeType>
void
ComputeFVFaceResidualsThread<RangeType>::pre()
{
}

template <typename RangeType>
void
ComputeFVFaceResidualsThread<RangeType>::post()
{
}

template <typename RangeType>
void
ComputeFVFaceResidualsThread<RangeType>::onBoundary(const FaceInfo & fi, BoundaryID bnd_id)
{
  std::vector<FVBoundaryCondition *> bcs;
  _fe_problem.theWarehouse()
      .query()
      .condition<AttribSystem>("FVBC")
      .condition<AttribBoundary>(bnd_id)
      .queryInto(bcs);
  if (bcs.size() == 0)
    return;

  _fe_problem.reinitElemFace(fi.leftElem(), fi.leftSide(), bnd_id, _tid);

  // Set up Sentinel class so that, even if reinitMaterialsFace() throws, we
  // still remember to swap back during stack unwinding.
  SwapBackSentinel sentinel(_fe_problem, &FEProblem::swapBackMaterialsFace, _tid);

  _fe_problem.reinitMaterialsFace(elem->subdomain_id(), _tid);
  _fe_problem.reinitMaterialsBoundary(bnd_id, _tid);

  for (const auto & bc : bcs)
    bc->computeResidual();
}

template <typename RangeType>
void
ComputeFVFaceResidualsThread<RangeType>::onFace(const FaceInfo & fi)
{
  std::vector<FVKernel *> kernels;
  _fe_problem.theWarehouse()
      .query()
      .condition<AttribSystem>("FVFluxKernels")
      .condition<AttribSubdomains>(_subdomain)
      .queryInto(kernels);
  if (kernels.size() == 0)
    return;

  // TODO: we should do reconstruction on-demand inside the FV
  // variable class and return reconstructed values for face reinits and
  // regular values for volume reinits.

  // This needs to skip all the regular FE and fancy quadrature (re)init stuff
  // and just update assembly and dof indices sorts of things as necessary.
  _fe_problem.reinitFVFace(fi);

  // Set up Sentinels so that, even if one of the reinitMaterialsXXX() calls throws, we
  // still remember to swap back during stack unwinding.
  SwapBackSentinel face_sentinel(_fe_problem, &FEProblem::swapBackMaterialsFace, _tid);
  _fe_problem.reinitMaterialsFace(fi.leftElem().subdomain_id(), _tid);

  // this reinits the materials for the neighbor element on the face (although it may not look like
  // it)
  SwapBackSentinel neighbor_sentinel(_fe_problem, &FEProblem::swapBackMaterialsNeighbor, _tid);
  _fe_problem.reinitMaterialsNeighbor(fi.rightElem().subdomain_id(), _tid);

  for (const auto k : kernels)
    auto r = k->computeResidual();
}

template <typename RangeType>
void
ComputeFVFaceResidualsThread<RangeType>::subdomainChanged()
{
  std::set<MooseVariableFEBase *> needed_mat_props;
  std::set<unsigned int> needed_mat_props;

  // TODO: do this for other relevant objects - like FV BCs, etc.
  std::vector<FVKernel *> kernels;
  _fe_problem.theWarehouse()
      .query()
      .condition<AttribSystem>("FVFluxKernels")
      .condition<AttribSubdomains>(_subdomain)
      .queryInto(kernels);
  for (auto k : kernels)
  {
    const auto & deps = k->getMooseVariableDependencies();
    needed_moose_vars.insert(deps.begin(), deps.end());
    const auto & mdeps = object->getMatPropDependencies();
    needed_mat_props.insert(mdeps.begin(), mdeps.end());
  }

  _fe_problem.setActiveElementalMooseVariables(needed_moose_vars, _tid);
  _fe_problem.setActiveMaterialProperties(needed_mat_props, _tid);
  _fe_problem.prepareMaterials(_subdomain, _tid);
}

template <typename RangeType>
void
ComputeFVFaceResidualsThread<RangeType>::neighborSubdomainChanged()
{
}
