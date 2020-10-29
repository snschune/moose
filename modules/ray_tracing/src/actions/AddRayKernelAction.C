//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AddRayKernelAction.h"

// Local includes
#include "RayKernelBase.h"

registerMooseAction("RayTracingApp", AddRayKernelAction, "add_ray_kernel");

InputParameters
AddRayKernelAction::validParams()
{
  return AddRayTracingObjectAction::validParams();
}

AddRayKernelAction::AddRayKernelAction(InputParameters params) : AddRayTracingObjectAction(params)
{
}

void
AddRayKernelAction::act()
{
  setRayTracingStudy();
  _problem->addObject<RayKernelBase>(_type, _name, _moose_object_pars);
}
