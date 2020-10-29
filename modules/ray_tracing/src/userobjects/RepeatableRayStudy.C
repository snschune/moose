//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RepeatableRayStudy.h"

registerMooseObject("RayTracingApp", RepeatableRayStudy);

InputParameters
RepeatableRayStudy::validParams()
{
  auto params = RepeatableRayStudyBase::validParams();

  params.addRequiredParam<std::vector<std::string>>("names", "Unique names for the Rays");

  params.addRequiredParam<std::vector<Point>>("start_points", "The points to start Rays from");
  params.addParam<std::vector<Point>>("directions",
                                      "The directions to spawn Rays in (they do not need to be "
                                      "normalized to 1). Use either this parameter "
                                      "or the end_points parameter, but not both!");
  params.addParam<std::vector<Point>>("end_points",
                                      "The points where Rays should end. Use either this parameter "
                                      "or the directions parameter, but not both!");
  params.addParam<std::vector<Real>>(
      "max_distances",
      "The maximum distances that each Ray can travel before it is killed internally after "
      "RayKernel execution. This can ONLY be used when the 'directions' parameter is used. "
      "When this not set, the Rays must be killed by RayKernels, RayBCs, or the global "
      " max distance parameter, 'ray_max_distance' (applies to all Rays)");

  params.addParamNamesToGroup("start_points directions end_points max_distances", "Trajectory");

  params.addParam<std::vector<std::string>>(
      "ray_data_names",
      "The Ray data names to register. If 'initial_ray_data' is set, these data names will be "
      "associated with said initial values. Otherwise, they will be set to zero.");
  params.addParam<std::vector<Real>>(
      "initial_ray_data",
      "The initial Ray data to set. This must be paired with 'ray_data_names' to know what data "
      "names to assign the data to");
  params.addParam<std::vector<std::string>>(
      "ray_aux_data_names",
      "The Ray aux data names to register. If 'initial_ray_aux_data' is set, these aux data names "
      "will be associated with said initial values. Otherwise, they will be set to zero.");
  params.addParam<std::vector<Real>>(
      "initial_ray_aux_data",
      "The initial Ray aux data to set. This must be paired with 'ray_aux_data_names' to know what "
      "aux data names to assign the data to");

  params.addParamNamesToGroup(
      "ray_data_names initial_ray_data ray_aux_data_names initial_ray_aux_data", "Data");

  return params;
}

RepeatableRayStudy::RepeatableRayStudy(const InputParameters & parameters)
  : RepeatableRayStudyBase(parameters),
    _ids(registerRays(getParam<std::vector<std::string>>("names"), /* graceful = */ true)),
    _start_points(getParam<std::vector<Point>>("start_points")),
    _end_points(parameters.isParamSetByUser("end_points")
                    ? &getParam<std::vector<Point>>("end_points")
                    : nullptr),
    _directions(parameters.isParamSetByUser("directions")
                    ? &getParam<std::vector<Point>>("directions")
                    : nullptr),
    _max_distances(parameters.isParamSetByUser("max_distances")
                       ? &getParam<std::vector<Real>>("max_distances")
                       : nullptr),
    _ray_data_indices(parameters.isParamSetByUser("ray_data_names")
                          ? registerRayData(getParam<std::vector<std::string>>("ray_data_names"))
                          : std::vector<RayDataIndex>()),
    _initial_ray_data(parameters.isParamSetByUser("initial_ray_data")
                          ? &getParam<std::vector<Real>>("initial_ray_data")
                          : nullptr),
    _ray_aux_data_indices(
        parameters.isParamSetByUser("ray_aux_data_names")
            ? registerRayAuxData(getParam<std::vector<std::string>>("ray_aux_data_names"))
            : std::vector<RayDataIndex>()),
    _initial_ray_aux_data(parameters.isParamSetByUser("initial_ray_aux_data")
                              ? &getParam<std::vector<Real>>("initial_ray_aux_data")
                              : nullptr)
{
  for (std::size_t i = 0; i < _ids.size(); ++i)
    if (_ids[i] == Ray::INVALID_RAY_ID)
      paramError("names",
                 "Ray name '",
                 getParam<std::string>("names")[i],
                 "' could not be registered. It may have been provided multiple times.");
  if (_end_points && _directions)
    paramError("directions", "Can only use 'directions' or 'end_points', but not both");
  if (_start_points.size() != _ids.size())
    paramError("start_points", "Not the same size as names");
  if (_directions && _ids.size() != _directions->size())
    paramError("directions", "Not the same size as names");
  if (_max_distances)
  {
    if (!_directions)
      paramError("max_distances",
                 "Can only be used when trajectories are set with the 'directions' parameter");
    if (_max_distances->size() != _start_points.size())
      paramError("max_distances", "Must be the same size as 'start_points'");
    for (const auto val : *_max_distances)
      if (val <= 0)
        paramError("max_distances", "Values must be positive");
  }
  if (_end_points && _ids.size() != _end_points->size())
    paramError("end_points", "Not the same size as names");
  if (!_ray_data_indices.empty() && _initial_ray_data &&
      _initial_ray_data->size() != _ray_data_indices.size())
    paramError("initial_ray_data", "Not the same size as ray_data_names");
  if (_initial_ray_data && _ray_data_indices.empty())
    paramError("initial_ray_data", "Can only be used if ray_data_names is set");
  if (!_ray_aux_data_indices.empty() && _initial_ray_aux_data &&
      _initial_ray_aux_data->size() != _ray_aux_data_indices.size())
    paramError("initial_ray_aux_data", "Not the same size as ray_aux_data_names");
  if (_initial_ray_aux_data && _ray_aux_data_indices.empty())
    paramError("initial_ray_aux_data", "Can only be used if ray_aux_data_names is set");
}

void
RepeatableRayStudy::defineRays()
{
  for (std::size_t i = 0; i < _ids.size(); ++i)
  {
    std::shared_ptr<Ray> ray = acquireRay(/* tid = */ 0);

    if (_end_points) // user set end point
      ray->setStartEnd(_start_points[i], (*_end_points)[i]);
    else // user set direction
      ray->setStartDirection(_start_points[i], (*_directions)[i]);

    // We already generated the IDs
    ray->setID(_ids[i]);

    // Set the data
    ray->data().resize(rayDataSize(), 0);
    if (_initial_ray_data)
      for (std::size_t i = 0; i < _ray_data_indices.size(); ++i)
        ray->data(_ray_data_indices[i]) = (*_initial_ray_data)[i];

    // Set the aux data
    ray->auxData().resize(rayAuxDataSize(), 0);
    if (_initial_ray_aux_data)
      for (std::size_t i = 0; i < _ray_aux_data_indices.size(); ++i)
        ray->auxData(_ray_aux_data_indices[i]) = (*_initial_ray_aux_data)[i];

    // User set max-distances
    if (_max_distances)
      ray->setMaxDistance((*_max_distances)[i]);

    _rays.emplace_back(std::move(ray));
  }
}
