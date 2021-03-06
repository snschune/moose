/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "ComputeElemAuxVarsThread.h"
#include "AuxiliarySystem.h"
#include "AuxKernel.h"
#include "FEProblem.h"
// libmesh includes
#include "libmesh/threads.h"


ComputeElemAuxVarsThread::ComputeElemAuxVarsThread(FEProblem & problem, AuxiliarySystem & sys, std::vector<AuxWarehouse> & auxs) :
    ThreadedElementLoop<ConstElemRange>(problem, sys),
    _aux_sys(sys),
    _auxs(auxs)
{
}

// Splitting Constructor
ComputeElemAuxVarsThread::ComputeElemAuxVarsThread(ComputeElemAuxVarsThread & x, Threads::split /*split*/) :
    ThreadedElementLoop<ConstElemRange>(x._fe_problem, x._system),
    _aux_sys(x._aux_sys),
    _auxs(x._auxs)
{
}

ComputeElemAuxVarsThread::~ComputeElemAuxVarsThread()
{
}

void
ComputeElemAuxVarsThread::subdomainChanged()
{
  // prepare variables
  for (std::map<std::string, MooseVariable *>::iterator it = _aux_sys._elem_vars[_tid].begin(); it != _aux_sys._elem_vars[_tid].end(); ++it)
  {
    MooseVariable * var = it->second;
    var->prepareAux();
  }

  // block setup
  for(std::vector<AuxKernel *>::const_iterator aux_it=_auxs[_tid].activeBlockElementKernels(_subdomain).begin();
      aux_it != _auxs[_tid].activeBlockElementKernels(_subdomain).end();
      aux_it++)
    (*aux_it)->subdomainSetup();

  std::set<MooseVariable *> needed_moose_vars;

  for(std::vector<AuxKernel*>::const_iterator block_element_aux_it = _auxs[_tid].activeBlockElementKernels(_subdomain).begin();
      block_element_aux_it != _auxs[_tid].activeBlockElementKernels(_subdomain).end(); ++block_element_aux_it)
  {
    const std::set<MooseVariable *> & mv_deps = (*block_element_aux_it)->getMooseVariableDependencies();
    needed_moose_vars.insert(mv_deps.begin(), mv_deps.end());
  }

  _fe_problem.setActiveElementalMooseVariables(needed_moose_vars, _tid);
  _fe_problem.prepareMaterials(_subdomain, _tid);
}


void
ComputeElemAuxVarsThread::onElement(const Elem * elem)
{
  if (! _auxs[_tid].activeBlockElementKernels(_subdomain).empty())
  {
    _fe_problem.prepare(elem, _tid);
    _fe_problem.reinitElem(elem, _tid);
    _fe_problem.reinitMaterials(elem->subdomain_id(), _tid);

    for(std::vector<AuxKernel*>::const_iterator block_element_aux_it = _auxs[_tid].activeBlockElementKernels(_subdomain).begin();
        block_element_aux_it != _auxs[_tid].activeBlockElementKernels(_subdomain).end(); ++block_element_aux_it)
      (*block_element_aux_it)->compute();

    _fe_problem.swapBackMaterials(_tid);

    // update the solution vector
    {
      Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
      for (std::map<std::string, MooseVariable *>::iterator it = _aux_sys._elem_vars[_tid].begin(); it != _aux_sys._elem_vars[_tid].end(); ++it)
      {
        MooseVariable * var = it->second;
        var->insert(_system.solution());
      }
    }
  }
}

void
ComputeElemAuxVarsThread::post()
{
  _fe_problem.clearActiveElementalMooseVariables(_tid);
}

void
ComputeElemAuxVarsThread::join(const ComputeElemAuxVarsThread & /*y*/)
{
}
