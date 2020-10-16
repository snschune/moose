//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AddLotsOfFVAdvectionReaction.h"
#include "Parser.h"
#include "FEProblem.h"
#include "Factory.h"
#include "MooseEnum.h"
#include "AddVariableAction.h"
#include "Conversion.h"
#include "DirichletBC.h"
#include "AddVariableAction.h"

#include <sstream>
#include <stdexcept>

#include "libmesh/libmesh.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/explicit_system.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/fe.h"

registerMooseAction("MooseTestApp", AddLotsOfFVAdvectionReaction, "add_variable");
registerMooseAction("MooseTestApp", AddLotsOfFVAdvectionReaction, "add_fv_kernel");
registerMooseAction("MooseTestApp", AddLotsOfFVAdvectionReaction, "add_fv_bc");
registerMooseAction("MooseTestApp", AddLotsOfFVAdvectionReaction, "add_ic");

InputParameters
AddLotsOfFVAdvectionReaction::validParams()
{
  InputParameters params = Action::validParams();
  params.addRequiredParam<RealVectorValue>("vel", "Advecting velocity");
  params.addRequiredParam<std::vector<BoundaryName>>("inlet_boundaries",
                                                     "The list of inlet boundary IDs");
  params.addRequiredParam<std::vector<BoundaryName>>("outlet_boundaries",
                                                     "The list of outlet boundary IDs");
  params.addRequiredParam<unsigned int>("number", "The number of variables to add");
  return params;
}

AddLotsOfFVAdvectionReaction::AddLotsOfFVAdvectionReaction(const InputParameters & params)
  : Action(params)
{
}

VariableName
AddLotsOfFVAdvectionReaction::varName(unsigned int j) const
{
  std::stringstream ss;
  ss << "advected_var_" << j;
  return ss.str();
}

void
AddLotsOfFVAdvectionReaction::act()
{
  unsigned int number = getParam<unsigned int>("number");

  if (_current_task == "add_variable")
  {
    // first tell the FEProblem that fv is required
    _problem->needFV();

    // define common variable parameters
    const FEType fe_type(Utility::string_to_enum<Order>("CONSTANT"),
                         Utility::string_to_enum<FEFamily>("MONOMIAL"));
    auto type = AddVariableAction::determineType(fe_type, 1, true);
    auto var_params = _factory.getValidParams(type);
    var_params.set<MooseEnum>("family") = Moose::stringify(fe_type.family);
    var_params.set<MooseEnum>("order") = fe_type.order.get_order();
    var_params.set<bool>("fv") = true;

    // add pebble speed variable
    for (unsigned int j = 0; j < number; ++j)
      _problem->addVariable(type, varName(j), var_params);
  }
  else if (_current_task == "add_ic")
  {
    for (unsigned int j = 0; j < number; ++j)
    {
      std::stringstream ss;
      ss << "ic_" << j;
      InputParameters params = _factory.getValidParams("ConstantIC");
      params.set<VariableName>("variable") = varName(j);
      params.set<Real>("value") = 1;
      _problem->addInitialCondition("ConstantIC", ss.str(), params);
    }
  }
  else if (_current_task == "add_fv_kernel")
  {
    for (unsigned int j = 0; j < number; ++j)
    {
      std::stringstream ss;
      ss << "advection_" << j;
      InputParameters params = _factory.getValidParams("FVAdvection");
      params.set<NonlinearVariableName>("variable") = varName(j);
      params.set<RealVectorValue>("velocity") = getParam<RealVectorValue>("vel");
      _problem->addFVKernel("FVAdvection", ss.str(), params);
    }

    // reaction term for all but the last
    for (unsigned int j = 0; j < number - 1; ++j)
    {
      std::stringstream ss;
      ss << "reaction_" << j;
      InputParameters params = _factory.getValidParams("FVReaction");
      params.set<NonlinearVariableName>("variable") = varName(j);
      _problem->addFVKernel("FVReaction", ss.str(), params);
    }

    // a coupled source for all but the first
    for (unsigned int j = 1; j < number; ++j)
    {
      std::stringstream ss;
      ss << "coupled_" << j;
      InputParameters params = _factory.getValidParams("FVCoupledForce");
      params.set<NonlinearVariableName>("variable") = varName(j);
      params.set<std::vector<VariableName>>("v") = {varName(j - 1)};
      _problem->addFVKernel("FVCoupledForce", ss.str(), params);
    }
  }
  else if (_current_task == "add_fv_bc")
  {
    // only the first variables gets a nonzero inflow
    InputParameters params_inlet = _factory.getValidParams("FVNeumannBC");
    params_inlet.set<NonlinearVariableName>("variable") = varName(0);
    params_inlet.set<Real>("value") = 1;
    params_inlet.set<std::vector<BoundaryName>>("boundary") =
        getParam<std::vector<BoundaryName>>("inlet_boundaries");
    _problem->addFVBC("FVNeumannBC", "inlet", params_inlet);

    // outlet
    for (unsigned int j = 0; j < number; ++j)
    {
      std::stringstream ss;
      ss << "outflow_" << j;
      InputParameters params = _factory.getValidParams("FVConstantScalarOutflowBC");
      params.set<NonlinearVariableName>("variable") = varName(j);
      params.set<RealVectorValue>("velocity") = getParam<RealVectorValue>("vel");
      params.set<std::vector<BoundaryName>>("boundary") =
          getParam<std::vector<BoundaryName>>("outlet_boundaries");
      _problem->addFVBC("FVConstantScalarOutflowBC", ss.str(), params);
    }
  }
}
