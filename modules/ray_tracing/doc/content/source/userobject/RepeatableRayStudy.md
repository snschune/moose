# RepeatableRayStudy

## Description

The `RepeatableRayStudy` is a specialized [RayTracingStudy.md] that generates rays from a set of user-input start points and end points/directions.

!alert note
The `RepeatableRayStudy` is meant to be the primary [RayTracingStudy.md] to be used for the majority of use cases and is not meant to be derived from! It is "repeatable" because it works with adaptivity and multiple executions (transients, on residual/Jacobian evaluations, etc). It does not require the user to have any knowledge of how to generate rays or determine on which processor and element element rays need to start depending on their starting point.

## Defining the Rays

The following parameters must be set:

- `names`: A list of unique names to identify the rays being generated.
- `start_points`: A list of points that the rays should start from.

When using [RayKernels/index.md] and [RayBCs/index.md] with the `RepeatableRayStudy`, you must specify which rays the [RayKernels/index.md]/[RayBCs/index.md] are applied to via their own `names` parameter. The names supplied to the [RayKernels/index.md] and [RayBCs/index.md] are the same as the names specified in your `RepeatableRayStudy`.

After setting these parameters, you must decide if you want to define the remainder of the trajectory by end points or by directions.

### Defining By End Points

To define the remainder of the trajectory by end points, provide the points at which you want the rays to end in the `end_points` parameter. When the [Ray.md] end points are set, internally the tracer will set the max distance of each [Ray.md] individually such that they all end at the straight-line distance between the provided start point and the provided end point.

[RayKernels/index.md] and [RayBCs/index.md] can still kill the rays earlier along their trajectory, but they are guaranteed to end once they hit either their end point or possibly sooner if the study's global maximum ray distance (the `ray_distance` parameter) is less than the distance from the start to end point.

Rays that are killed due to reaching their max distance (which is the case when they reach their end point) are killed before the execution of [RayBCs/index.md]. For example, if a [Ray.md] reaches its end point and said end point is on a boundary with [RayBCs/index.md], the [RayBCs/index.md] will not be executed on the [Ray.md].

!alert note
Rays that have had their trajectory set via end points are not allowed to have their trajectories modified mid-trace via [RayKernels/index.md] or [RayBCs/index.md]. For example, these rays cannot be reflected on boundaries via the [ReflectRayBC.md]. You must instead define rays by the `directions` parameter if you want them to be able to have their trajectories changed mid-trace.

### Defining By Directions

To define the remainder of the trajectory by directions, provide the directions at which you want the rays to travel in the `directions` parameter. These directions do not need to be normalized.

When the [Ray.md] trajectory is defined by a direction, the user is responsible for killing the [Ray.md]. If the [Ray.md] hits an external boundary and has not been killed or had its trajectory changed, it will error.

Common ways of ending a [Ray.md] when its trajectory is defined by a direction:

- Killing the [Ray.md] on a boundary via [RayBCs/index.md] (see [KillRayBC.md]).
- Setting the maximum distance each [Ray.md] can travel via the `max_distances` parameter.
- Setting the maximum distance all [Ray.md]s can travel via the `ray_distance` parameter.
- Killing the [Ray.md] within an element via [RayKernels/index.md] (see [RayKernels/index.md#ending-the-ray]).

## Setting Ray Data

For more advanced use, one can also register [Ray.md] data/auxiliary data and initialize it as desired. It is important that this is not necessary when using [RayKernels/index.md] that contribute to residuals or integrate along lones, as the [Ray.md] data mangement in those cases is handled under the hood.

!syntax parameters /UserObjects/RepeatableRayStudy
  visible=Required Trajectory

!syntax inputs /UserObjects/RepeatableRayStudy
