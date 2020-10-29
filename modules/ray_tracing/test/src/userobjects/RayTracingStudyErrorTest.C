//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RayTracingStudyErrorTest.h"

registerMooseObject("RayTracingTestApp", RayTracingStudyErrorTest);

InputParameters
RayTracingStudyErrorTest::validParams()
{
  auto params = RayTracingStudy::validParams();

  params.addParam<bool>(
      "add_outside_of_generation", false, "True to add a Ray outside of generation");
  params.addParam<bool>(
      "add_thread1_during_generation", false, "True to add a Ray during generation on thread 1");
  params.set<bool>("_use_ray_registration") = false;

  return params;
}

RayTracingStudyErrorTest::RayTracingStudyErrorTest(const InputParameters & parameters)
  : RayTracingStudy(parameters)
{
}

void
RayTracingStudyErrorTest::initialSetup()
{
  RayTracingStudy::initialSetup();

  auto ray = acquireRay(0);
  ray->setID(generateUniqueRayID(0));
  ray->setStartingElem(_mesh.elemPtr(0));
  ray->setStartDirection(_mesh.elemPtr(0)->centroid(), Point(1, 1, 0));

  if (getParam<bool>("add_outside_of_generation"))
    moveRayToBuffer(ray);
}

void
RayTracingStudyErrorTest::generateRays()
{
  auto ray = acquireRay(0);
  ray->setID(generateUniqueRayID(0));
  ray->setStartingElem(_mesh.elemPtr(0));
  ray->setStartDirection(_mesh.elemPtr(0)->centroid(), Point(1, 1, 0));

  if (getParam<bool>("add_thread1_during_generation"))
    moveRayToBuffer(ray, 1);
}
