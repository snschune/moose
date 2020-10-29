#pragma once

#include "ParallelStudy.h"

// Local includes
#include "Ray.h"

// Forward declarations
class RayTracingStudy;
class TraceRay;

class ParallelRayStudy : public ParallelStudy<Ray, RayTracingStudy>
{
public:
  ParallelRayStudy(RayTracingStudy & study,
                   const std::vector<std::shared_ptr<TraceRay>> & threaded_trace_ray);

protected:
  void executeObject(const std::shared_ptr<Ray> & ray, const THREAD_ID tid) override;
  bool shouldExecuteObject(const std::shared_ptr<Ray> & object) const override;
  bool shouldCommunicateObject(const std::shared_ptr<Ray> & ray,
                               processor_id_type & to_pid) const override;
  void onCompleteObject(const std::shared_ptr<Ray> & ray) override;
  std::string addObjectError(const std::shared_ptr<Ray> & ray,
                             const AddObjectError error) const override;

  /// The RayTracingStudy
  RayTracingStudy & _ray_tracing_study;
  /// The TraceRay objects that do the tracing for each thread
  const std::vector<std::shared_ptr<TraceRay>> & _threaded_trace_ray;
};
