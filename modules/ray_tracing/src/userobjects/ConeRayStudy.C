//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ConeRayStudy.h"

// Local includes
#include "RayTracingAngularQuadrature.h"

registerMooseObject("RayTracingApp", ConeRayStudy);

InputParameters
ConeRayStudy::validParams()
{
  auto params = RepeatableRayStudyBase::validParams();

  params.addClassDescription("Ray study that spawns Rays in a cone from a given set of starting "
                             "points for the cones and half angles for the cones.");

  params.addRequiredParam<std::vector<Point>>("start_points", "The point(s) of the cone(s).");
  params.addRequiredParam<std::vector<Point>>(
      "directions", "The direcions(s) of the cone(s) (points down the center of each cone).");
  params.addParam<std::vector<Real>>("scaling_factors",
                                     "Scaling factors for each cone (if any). Defaults to 1.");

  params.addRequiredParam<std::vector<Real>>("half_cone_angles",
                                             "Angle of the half-cones in degrees (must be <= 90)");
  params.addParam<std::vector<unsigned int>>(
      "polar_quad_orders",
      "Order of the polar quadrature for each cone. Polar angle is between ray and the "
      "cone direction. Must be odd. This will default to 15 for all cones if not provided.");
  params.addParam<unsigned int>(
      "azimuthal_quad_orders",
      "Order of the azimuthal quadrature per quadrant for each cone. The azimuthal angle is "
      "measured in a plane perpendicular to the cone direction. Not needed in 2D. This will "
      "default to 8 if not provided.");

  params.addRequiredParam<std::string>("ray_data_name",
                                       "The name of the Ray data that the angular quadrature "
                                       "weights and factors will be filled into for "
                                       "properly weighting the line source per Ray.");

  // Here we will not use Ray registration because we create so many Rays that
  // we don't want to keep track of them by name
  params.set<bool>("_use_ray_registration") = false;

  return params;
}

ConeRayStudy::ConeRayStudy(const InputParameters & parameters)
  : RepeatableRayStudyBase(parameters),
    _start_points(getParam<std::vector<Point>>("start_points")),
    _directions(getParam<std::vector<Point>>("directions")),
    _scaling_factors(parameters.isParamSetByUser("scaling_factors")
                         ? getParam<std::vector<Real>>("scaling_factors")
                         : std::vector<Real>(_start_points.size(), 1)), // default to 1
    _half_cone_angles(getParam<std::vector<Real>>("half_cone_angles")),
    _polar_quad_orders(parameters.isParamSetByUser("polar_quad_orders")
                           ? getParam<std::vector<unsigned int>>("polar_quad_orders")
                           : std::vector<unsigned int>(_start_points.size(), 15)), // default to 15
    _azimuthal_quad_orders(
        parameters.isParamSetByUser("azimuthal_quad_orders")
            ? getParam<std::vector<unsigned int>>("azimuthal_quad_orders")
            : std::vector<unsigned int>(_start_points.size(), 8)), // default to 8
    _ray_data_index(registerRayData(getParam<std::string>("ray_data_name")))
{
  if (_directions.size() != _start_points.size())
    paramError("directions", "Not the same size as start_points.");

  if (_scaling_factors.size() != _start_points.size())
    paramError("scaling_factors", "Not the same size as start_points.");

  if (_half_cone_angles.size() != _start_points.size())
    paramError("half_cone_angles", "Not the same size as start_points.");
  for (const auto val : _half_cone_angles)
    if (val <= 0 || val > 90)
      paramError("half_cone_angles", "Must be > 0 and <= 90 degrees");

  if (_polar_quad_orders.size() != _start_points.size())
    paramError("polar_quad_orders", "Not the same size as start_points.");
  for (const auto val : _polar_quad_orders)
    if (val % 2 == 0)
      paramError("polar_quad_orders", "Must be odd.");

  if (_azimuthal_quad_orders.size() != _start_points.size())
    paramError("azimuthal_quad_orders", "Not the same size as start_points.");
  if (_mesh.dimension() == 2 && parameters.isParamSetByUser("azimuthal_quad_orders"))
    paramError("azimuthal_quad_orders", "Not required for 2D.");
  for (const auto val : _azimuthal_quad_orders)
    if (val == 0)
      paramError("azimuthal_quad_orders", "Must be nonzero.");

  if (_mesh.dimension() == 1)
    mooseError(_error_prefix, ": Does not support 1D.");
}

void
ConeRayStudy::defineRays()
{
  // Loop through each cone
  for (std::size_t i = 0; i < _start_points.size(); ++i)
  {
    // The half cone angle in rad
    const auto half_cone_rad = _half_cone_angles[i] * M_PI / 180.0;

    // Create the angular quadrature for this cone
    std::vector<std::pair<Real, Real>> angles;
    std::vector<Real> weights;
    if (_mesh.dimension() == 2)
      RayTracingAngularQuadrature::getHalfRange2D(
          _polar_quad_orders[i], angles, weights, half_cone_rad);
    else
      RayTracingAngularQuadrature::getHalfRange3D(
          4 * _azimuthal_quad_orders[i], _polar_quad_orders[i], angles, weights, half_cone_rad);

    // This rotation matrix is needed to rotate the angular quadrature about
    // the reference direction (the direction that points down the center of the cone)
    const auto rotation_matrix =
        RayTracingAngularQuadrature::getRotationMatrix(_directions[i], _mesh.dimension());

    // For all angles in the angular quadrature, spawn a Ray
    for (std::size_t l = 0; l < angles.size(); ++l)
    {
      // The angular quadrature direction rotated about the reference direction
      // that is down the middle of the cone
      const auto direction =
          RayTracingAngularQuadrature::getDirection(l, rotation_matrix, angles, _mesh.dimension());

      // Get a Ray from the study to initialize
      std::shared_ptr<Ray> ray = acquireRay(/* tid = */ 0);

      // Start from our cone point in the rotated angular quadrature direction
      // Note here that we do not need to set the starting element - all Rays
      // at this point are replicated across all processors and will be
      // properly claimed (moved to the starting proc with the correct starting elem)
      ray->setStartDirection(_start_points[i], direction);

      // Unique ID (but replicated across all processors) for this point/direction combo
      ray->setID(i * angles.size() + l);

      // Size and zero all of the data
      ray->data().resize(rayDataSize(), 0);
      ray->auxData().resize(rayAuxDataSize(), 0);

      // Add the angular quadrature weight and scaling factor as data on the Ray for weighting
      // The weights currently sum to 2 * half_cone_rad - we want them to sum to 1
      ray->data(_ray_data_index) = _scaling_factors[i] * weights[l] / (2 * half_cone_rad);

      // Done with this Ray - move it to be traced later on
      _rays.emplace_back(std::move(ray));
    }
  }
}
