[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 5
    ny = 5
    xmax = 5
    ymax = 5
  []
[]

[Variables/u]
[]

[Kernels/diff]
  type = Diffusion
  variable = u
[]

[BCs]
  [left]
    type = DirichletBC
    variable = u
    boundary = left
    value = 0
  []
  [right]
    type = DirichletBC
    variable = u
    boundary = right
    value = 1
  []
[]

[Executioner]
  type = Steady
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  exodus = false
  csv = true
[]

[UserObjects/study]
  type = RepeatableRayStudy
  names = 'diag
           right_up'
  start_points = '0 0 0
                  5 0 0'
  end_points = '5 5 0
                5 5 0'
  execute_on = TIMESTEP_END
[]

[RayKernels/u_integral]
  type = VariableIntegralRayKernel
  variable = u
  rays = 'diag right_up'
[]

[Postprocessors]
  [diag_line_integral]
    type = RayIntegralValue
    ray_kernel = u_integral
    ray = diag
  []
  [right_up_line_integral]
    type = RayIntegralValue
    ray_kernel = u_integral
    ray = right_up
  []
[]
