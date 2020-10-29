[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 10
    ny = 10
    xmax = 5
    ymax = 5
  []
[]

[Variables/u]
[]

[Kernels]
  [reaction]
    type = Reaction
    variable = u
  []
  [diffusion]
    type = Diffusion
    variable = u
  []
[]

[UserObjects/study]
  type = ConeRayStudy
  start_points = '1 1.5 0'
  directions = '2 1 0'
  half_cone_angles = 2

  # Must be set with RayKernels that
  # contribute to the residual
  execute_on = PRE_KERNELS

  # For outputting Rays
  # always_cache_traces = true

  ray_data_name = weight
[]

[RayKernels/null]
  type = NullRayKernel
[]

[RayBCs]
  [reflect]
    type = ReflectRayBC
    boundary = 'right'
  []
  [kill_rest]
    type = KillRayBC
    boundary = 'top left'
  []
  # RayBCs not needed on bottom
  # because Rays will never hit the bottom
[]

[RayKernels/line_source]
  type = LineSourceRayKernel
  variable = u

  # Scale by the weights in the ConeRayStudy
  ray_data_factor_names = weight
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  exodus = true

  # For outputting the Rays
  # To enable, set execute_on = FINAL
  [rays]
    type = RayTracingExodus
    study = study
    execute_on = NONE # FINAL
  []
[]

[Adaptivity]
  steps = 0 # 6 for pretty pictures
  marker = marker
  initial_marker = marker
  max_h_level = 6
  [Indicators/indicator]
    type = GradientJumpIndicator
    variable = u
  []
  [Markers/marker]
    type = ErrorFractionMarker
    indicator = indicator
    coarsen = 0.25
    refine = 0.5
  []
[]
