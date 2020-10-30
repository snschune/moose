#include "ParallelRayStudy.h"

// Local includes
#include "TraceRay.h"
#include "RayTracingStudy.h"

ParallelRayStudy::ParallelRayStudy(
    RayTracingStudy & ray_tracing_study,
    const std::vector<std::shared_ptr<TraceRay>> & threaded_trace_ray)
  : ParallelStudy<Ray, RayTracingStudy>(ray_tracing_study.comm(),
                                        ray_tracing_study,
                                        ray_tracing_study.parameters(),
                                        "ParallelRayStudy"),
    _ray_tracing_study(ray_tracing_study),
    _threaded_trace_ray(threaded_trace_ray)
{
}

bool
ParallelRayStudy::shouldExecuteObject(const std::shared_ptr<Ray> & ray) const
{
  // If this is false it means we have a banked Ray that is either supposed to go to another
  // processor or ended on the boundary for another processor
  return ray->currentElem()->processor_id() == _pid;
}

bool
ParallelRayStudy::shouldCommunicateObject(const std::shared_ptr<Ray> & ray,
                                          processor_id_type & to_pid) const
{
  // Communicate the Ray if it should continue and is going to another processor
  if (ray->shouldContinue())
  {
    const auto elem_pid = ray->currentElem()->processor_id();

    if (elem_pid != _pid)
    {
      to_pid = elem_pid;
      return true;
    }

    mooseError("Continuing Ray not going to another processor after trace\n\n",
               ray->getInfo(&_ray_tracing_study));
  }

  // Ray shouldn't continue so it shouldn't be communicated
  return false;
}

void
ParallelRayStudy::executeObject(const std::shared_ptr<Ray> & ray, const THREAD_ID tid)
{
  if (!ray->shouldContinue())
    mooseError("Tracing Ray that should not continue\n\n", ray->getInfo(&_ray_tracing_study));

  _threaded_trace_ray[tid]->trace(ray);
}

void
ParallelRayStudy::onCompleteObject(const std::shared_ptr<Ray> & ray)
{
  _ray_tracing_study.onCompleteRay(ray);
}

std::string
ParallelRayStudy::addObjectError(const std::shared_ptr<Ray> & ray, const AddObjectError error) const
{
  std::stringstream oss;
  oss << "In method " << _ray_tracing_study.type() << "::addRay(s)ToBuffer:\n";

  if (error == AddObjectError::ADDITION_DURING_EXECUTION_DISABLED)
  {
    oss << "Rays are being added to the buffer during propagation.\n\n";
    oss << "This is currently disabled.\n";
    oss << "To enable, set the parameter 'allow_addition_during_execution' to true.";
  }
  else if (error == AddObjectError::PRE_EXECUTION_AND_EXECUTION_ONLY)
    oss << "Rays can only be added to the buffer during generateRays() and tracing.";
  else if (error == AddObjectError::PRE_EXECUTION_ONLY)
    oss << "Rays can only be added to the buffer during generateRays().";
  else if (error == AddObjectError::PRE_EXECUTION_THREAD_0_ONLY)
    oss << "Rays can only be added on thread 0 during generateRays() (not thread safe)";
  else
    mooseError("Unknown AddObjectError");

  if (ray)
    oss << "\n\n" << ray->getInfo(&_ray_tracing_study);

  return oss.str();
}
