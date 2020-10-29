[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 1
    ny = 1
    xmax = 1
    ymax = 1
  []
[]

[UserObjects/study]
  type = RayTracingStudyErrorTest
[]

[RayKernels/null]
  type = NullRayKernel
[]

[Problem]
  solve = false
[]

[Executioner]
  type = Steady
[]
