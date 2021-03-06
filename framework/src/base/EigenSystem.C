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
#include "EigenSystem.h"
#include "KernelWarehouse.h"

#include "MaterialData.h"
#include "Factory.h"
#include "EigenKernel.h"

EigenSystem::EigenSystem(FEProblem & fe_problem, const std::string & name) :
    NonlinearSystem(fe_problem, name),
    _all_eigen_vars(false),
    _active_on_old(false),
    _sys_sol_old(NULL),
    _sys_sol_older(NULL),
    _aux_sol_old(NULL),
    _aux_sol_older(NULL)
{
}

EigenSystem::~EigenSystem()
{
}

void
EigenSystem::addKernel(const std::string & kernel_name, const std::string & name, InputParameters parameters)
{
  for (THREAD_ID tid = 0; tid < libMesh::n_threads(); tid++)
  {
    // Set the parameters for thread ID and material data
    parameters.set<THREAD_ID>("_tid") = tid;
    parameters.set<MaterialData *>("_material_data") = _fe_problem._material_data[tid];

    // In the case of EigenKernels, we might need to add two to the system
    if (parameters.have_parameter<bool>("eigen"))
    {
      {
        // EigenKernel
        parameters.set<bool>("implicit") = true;
        EigenKernel *ekernel = static_cast<EigenKernel *>(_factory.create(kernel_name, name, parameters));
        mooseAssert(ekernel != NULL, "Not an EigenKernel object");
        if (parameters.get<bool>("eigen")) markEigenVariable(parameters.get<NonlinearVariableName>("variable"));
        // Extract the SubdomainIDs from the object (via BlockRestrictable class)
        std::set<SubdomainID> blk_ids = ekernel->blockIDs();
        _kernels[tid].addKernel(ekernel, blk_ids);
        _fe_problem._objects_by_name[tid][name].push_back(ekernel);
      }
      if (parameters.get<bool>("eigen"))
      {
        // EigenKernel_old
        parameters.set<bool>("implicit") = false;
        std::string old_name(name + "_old");

        EigenKernel *ekernel = static_cast<EigenKernel *>(_factory.create(kernel_name, old_name, parameters));
        _eigen_var_names.insert(parameters.get<NonlinearVariableName>("variable"));
        // Extract the SubdomainIDs from the object (via BlockRestrictable class)
        std::set<SubdomainID> blk_ids = ekernel->blockIDs();
        _kernels[tid].addKernel(ekernel, blk_ids);
        _fe_problem._objects_by_name[tid][old_name].push_back(ekernel);
      }
    }
    else // Standard nonlinear system kernel
    {
      // Create the kernel object via the factory
      KernelBase *kernel = static_cast<KernelBase *>(_factory.create(kernel_name, name, parameters));
      mooseAssert(kernel != NULL, "Not a Kernel object");
      // Extract the SubdomainIDs from the object (via BlockRestrictable class)
      std::set<SubdomainID> blk_ids = kernel->blockIDs();
      _kernels[tid].addKernel(kernel, blk_ids);
      _fe_problem._objects_by_name[tid][name].push_back(kernel);
    }
  }
}

void
EigenSystem::markEigenVariable(const VariableName & var_name)
{
  _eigen_var_names.insert(var_name);
}

void
EigenSystem::scaleSystemSolution(SYSTEMTAG tag, Real scaling_factor)
{
  if (tag==ALL)
  {
    solution().scale(scaling_factor);
  }
  else if (tag==EIGEN)
  {
    if (_all_eigen_vars)
    {
      solution().scale(scaling_factor);
    }
    else
    {
      std::set<dof_id_type>::iterator it      = _eigen_var_indices.begin();
      std::set<dof_id_type>::iterator it_end  = _eigen_var_indices.end();

      for(; it !=it_end; ++it)
        solution().set( *it, solution()(*it)*scaling_factor );
    }
  }
  solution().close();
  update();
}

void
EigenSystem::combineSystemSolution(SYSTEMTAG tag, const std::vector<Real> & coefficients)
{
  mooseAssert(coefficients.size()>0 && coefficients.size()<=3, "Size error on coefficients");
  if (tag==ALL)
  {
    solution().scale(coefficients[0]);
    if (coefficients.size()>1) solution().add(coefficients[1], solutionOld());
    if (coefficients.size()>2) solution().add(coefficients[2], solutionOlder());
  }
  else if (tag==EIGEN)
  {
    if (_all_eigen_vars)
    {
      solution().scale(coefficients[0]);
      if (coefficients.size()>1) solution().add(coefficients[1], solutionOld());
      if (coefficients.size()>2) solution().add(coefficients[2], solutionOlder());
    }
    else
    {
      std::set<dof_id_type>::iterator it      = _eigen_var_indices.begin();
      std::set<dof_id_type>::iterator it_end  = _eigen_var_indices.end();

      if (coefficients.size()>2)
      {
        for(; it !=it_end; ++it)
        {
          Real t = solution()(*it) * coefficients[0];
          t += solutionOld()(*it) * coefficients[1];
          t += solutionOlder()(*it) * coefficients[2];
          solution().set( *it, t );
        }
      }
      else if (coefficients.size()>1)
      {
        for(; it !=it_end; ++it)
        {
          Real t = solution()(*it) * coefficients[0];
          t += solutionOld()(*it) * coefficients[1];
          solution().set( *it, t );
        }
      }
      else
      {
        for(; it !=it_end; ++it)
        {
          Real t = solution()(*it) * coefficients[0];
          solution().set( *it, t );
        }
      }
    }
  }
  solution().close();
  update();
}

void
EigenSystem::initSystemSolution(SYSTEMTAG tag, Real v)
{
  if (tag==ALL)
  {
    solution() = v;
  }
  else if (tag==EIGEN)
  {
    if (_all_eigen_vars)
    {
      solution() = v;
    }
    else
    {
      std::set<dof_id_type>::iterator it      = _eigen_var_indices.begin();
      std::set<dof_id_type>::iterator it_end  = _eigen_var_indices.end();

      for(; it !=it_end; ++it)
        solution().set( *it, v );
    }
  }
  solution().close();
  update();
}

void
EigenSystem::initSystemSolutionOld(SYSTEMTAG tag, Real v)
{
  if (tag==ALL)
  {
    solutionOld() = v;
  }
  else if (tag==EIGEN)
  {
    if (_all_eigen_vars)
    {
      solutionOld() = v;
    }
    else
    {
      std::set<dof_id_type>::iterator it      = _eigen_var_indices.begin();
      std::set<dof_id_type>::iterator it_end  = _eigen_var_indices.end();

      for(; it !=it_end; ++it)
        solutionOld().set( *it, v );
    }
  }
  solutionOld().close();
  update();
}

void
EigenSystem::eigenKernelOnOld()
{
  _active_on_old = true;
}

void
EigenSystem::eigenKernelOnCurrent()
{
  _active_on_old = false;
}

bool
EigenSystem::activeOnOld()
{
  return _active_on_old;
}

void
EigenSystem::buildSystemDoFIndices(SYSTEMTAG tag)
{
  if (tag==ALL)
  {
  }
  else if (tag==EIGEN)
  {
    // build DoF indices for the eigen system
    _eigen_var_indices.clear();
    _all_eigen_vars = getEigenVariableNames().size()==getVariableNames().size();
    if (!_all_eigen_vars)
    {
      for (std::set<VariableName>::const_iterator it=getEigenVariableNames().begin();
           it!=getEigenVariableNames().end(); it++)
      {
        unsigned int i = sys().variable_number(*it);
        sys().local_dof_indices(i, _eigen_var_indices);
      }
    }
  }
}

void
EigenSystem::saveOldSolutions()
{
  if (!_sys_sol_old)
    _sys_sol_old= &addVector("save_flux_old", false, PARALLEL);
  if (!_aux_sol_old)
    _aux_sol_old = &_fe_problem.getAuxiliarySystem().addVector("save_aux_old",  false, PARALLEL);
  if (!_sys_sol_older)
    _sys_sol_older = &addVector("save_flux_older", false, PARALLEL);
  if (!_aux_sol_older)
    _aux_sol_older = &_fe_problem.getAuxiliarySystem().addVector("save_aux_older",  false, PARALLEL);
  *_sys_sol_old   = solutionOld();
  *_sys_sol_older = solutionOlder();
  *_aux_sol_old   = _fe_problem.getAuxiliarySystem().solutionOld();
  *_aux_sol_older = _fe_problem.getAuxiliarySystem().solutionOlder();
}

void
EigenSystem::restoreOldSolutions()
{
  solutionOld() = *_sys_sol_old;
  solutionOlder() = *_sys_sol_older;
  _fe_problem.getAuxiliarySystem().solutionOld() = *_aux_sol_old;
  _fe_problem.getAuxiliarySystem().solutionOlder() = *_aux_sol_older;
}
