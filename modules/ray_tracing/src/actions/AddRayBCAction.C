//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AddRayBCAction.h"

// Local includes
#include "RayBoundaryConditionBase.h"

registerMooseAction("RayTracingApp", AddRayBCAction, "add_ray_boundary_condition");

InputParameters
AddRayBCAction::validParams()
{
  return AddRayTracingObjectAction::validParams();
}

AddRayBCAction::AddRayBCAction(InputParameters params) : AddRayTracingObjectAction(params) {}

void
AddRayBCAction::act()
{
  setRayTracingStudy();
  _problem->addObject<RayBoundaryConditionBase>(_type, _name, _moose_object_pars);
}
