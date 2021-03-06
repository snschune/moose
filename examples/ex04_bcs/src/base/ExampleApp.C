#include "ExampleApp.h"
#include "Moose.h"
#include "Factory.h"
#include "AppFactory.h"

// Example 4 Includes
#include "Convection.h"
#include "GaussContForcing.h"
#include "CoupledDirichletBC.h"
#include "CoupledNeumannBC.h"

template<>
InputParameters validParams<ExampleApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

ExampleApp::ExampleApp(const std::string & name, InputParameters parameters) :
    MooseApp(name, parameters)
{
  srand(processor_id());

  Moose::registerObjects(_factory);
  ExampleApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ExampleApp::associateSyntax(_syntax, _action_factory);
}

ExampleApp::~ExampleApp()
{
}

void
ExampleApp::registerObjects(Factory & factory)
{
  registerKernel(Convection);
  registerKernel(GaussContForcing);                 // Extra forcing term
  registerBoundaryCondition(CoupledDirichletBC);    // Register our Boundary Conditions
  registerBoundaryCondition(CoupledNeumannBC);
}

void
ExampleApp::registerApps()
{
  registerApp(ExampleApp);
}

void
ExampleApp::associateSyntax(Syntax& /*syntax*/, ActionFactory & /*action_factory*/)
{
}
