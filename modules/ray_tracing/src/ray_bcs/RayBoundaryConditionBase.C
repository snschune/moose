//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RayBoundaryConditionBase.h"

// Local includes
#include "RayTracingStudy.h"
#include "TraceRay.h"

// MOOSE includes
#include "Assembly.h"

InputParameters
RayBoundaryConditionBase::validParams()
{
  auto params = RayTracingObject::validParams();
  params += BoundaryRestrictableRequired::validParams();

  params.addParam<std::vector<std::string>>("depends_on",
                                            "Other RayBCs that this RayBC depends on");

  params.registerBase("RayBoundaryCondition");
  params.registerSystemAttributeName("RayBoundaryCondition");

  // We don't currently allow reinits on RayBCs just yet
  params.suppressParameter<bool>("implicit");

  return params;
}

RayBoundaryConditionBase::RayBoundaryConditionBase(const InputParameters & params)
  : RayTracingObject(params),
    Restartable(this, "RayBoundaryConditions"),
    BoundaryRestrictableRequired(this, false), // false for sidesets
    _current_intersection_point(_study.traceRay(_tid).currentIntersectionPoint()),
    _current_bnd_id(_study.traceRay(_tid).currentBoundaryID())
{
  // Add dependencies
  if (params.isParamSetByUser("depends_on"))
    for (const auto & name : getParam<std::vector<std::string>>("depends_on"))
      dependsOn(name);
}

RayBoundaryConditionBase::~RayBoundaryConditionBase() {}

void
RayBoundaryConditionBase::changeRayDirection(const Point & direction, const bool skip_changed_check)
{
  auto & ray = currentRay();

  if (!ray->shouldContinue())
    mooseError(_error_prefix,
               ": Cannot changeRayDirection() for a Ray that should not continue\n\n",
               ray->getInfo(&_study));

  if (!ray->intersections())
    mooseError(_error_prefix,
               ": Cannot change direction for Ray that has not moved\n\n",
               ray->getInfo(&_study));

  if (!skip_changed_check && ray->trajectoryChanged())
    mooseError(
        _error_prefix,
        " is trying to change a Ray's direction, but its direction has already been changed\n\n",
        ray->getInfo(&_study));

  if (ray->endSet())
    mooseError(_error_prefix,
               " is trying to change the direction of a Ray that"
               "\nhad its end point set upon generation (via setStartEnd())."
               "\n\nThis is not currently supported.\n\n",
               ray->getInfo(&_study));

  ray->changeDirection(direction, {});
}

std::shared_ptr<Ray>
RayBoundaryConditionBase::acquireRay(const Point & direction)
{
  mooseAssert(_study.currentlyPropagating(), "Should not be getting a Ray while not propagating");

  // Acquire a Ray with the proper sizes
  std::shared_ptr<Ray> ray = _study.acquireRay(_tid, _study.rayDataSize(), _study.rayAuxDataSize());
  // Set the starting element to the one that we're currently in
  ray->setStartingElem(_current_elem);
  // Set the incoming side to the side that we're currently on
  ray->setStartingIncomingSide(_current_intersected_side);
  // Set the start and direction with the current intersection and user specified direction
  ray->setStartDirection(_current_intersection_point, direction);
  // And give it a reasonable ID
  ray->setID(_study.generateUniqueRayID(_tid));

  return ray;
}

void
RayBoundaryConditionBase::moveRayToBuffer(std::shared_ptr<Ray> & ray)
{
  mooseAssert(_study.currentlyPropagating(),
              "Should not move Rays into buffer while not propagating");

  if (_current_elem != ray->startingElem())
    mooseError(_error_prefix,
               ": A Ray was added to the buffer mid-trace that does not\n",
               "start in the same Elem as the ",
               type(),
               " that created it\n\n",
               currentRay()->getInfo(&_study));

  if (ray->startingIncomingSide() != _current_intersected_side)
    mooseError(
        _error_prefix,
        " is trying to add a Ray to the buffer during tracing,\n",
        "but the incoming side of the new Ray is not the same as the intersected side of the ",
        type(),
        " that created it\n\n",
        ray->getInfo(&_study));

  _study.moveRayToBuffer(ray, _tid, _move_ray_error_prefix);
}
