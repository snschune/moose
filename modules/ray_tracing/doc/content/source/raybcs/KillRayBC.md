# KillRayBC

## Description

`KillRayBC` is a [RayBCs/index.md] that sets a [Ray.md] that is being traced to be killed on a boundary.

This is accomplished by setting the Ray to not be continued:

!listing modules/ray_tracing/src/ray_bcs/KillRayBC.C re=void\sKillRayBC::onBoundary.*?^}

!syntax parameters /RayBCs/KillRayBC

!syntax inputs /RayBCs/KillRayBC
