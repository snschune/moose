# RayDistanceAux

## Description

`RayDistanceAux` accumulates the distance traversed by each [Ray.md] segment into an `AuxVariable` for the element that the segment is in.

This is achieved by overriding `onSegment` and appending into the `AuxVariable` via `addValue()`, as:

!listing modules/ray_tracing/src/ray_kernels/RayDistanceAux.C re=void\sRayDistanceAux::onSegment.*?^}

!syntax parameters /RayKernels/RayDistanceAux

!syntax inputs /RayKernels/RayDistanceAux
