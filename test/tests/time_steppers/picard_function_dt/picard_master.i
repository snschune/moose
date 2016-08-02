[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 2
  ny = 2
  parallel_type = replicated
[]

[Variables]
  [./u]
  [../]
[]

[AuxVariables]
  [./v]
  [../]
[]

[Kernels]
  [./diff]
    type = CoefDiffusion
    variable = u
    coef = 0.1
  [../]
  [./time]
    type = TimeDerivative
    variable = u
  [../]
  [./force_u]
    type = CoupledForce
    variable = u
    v = v
  [../]
[]

[BCs]
  [./left]
    type = DirichletBC
    variable = u
    boundary = left
    value = 0
  [../]
  [./right]
    type = DirichletBC
    variable = u
    boundary = right
    value = 1
  [../]
[]

[Postprocessors]
  [./picard_its]
    type = NumPicardIterations
    execute_on = 'initial timestep_end'
  [../]
[]

[Executioner]
  # Preconditioned JFNK (default)
  type = Transient
  start_time = 0
  end_time = 0.9
  [./TimeStepper]
    type = PicardFunctionDT
    time_t  = '0 100'
    time_dt = '0.15 0.15'
    min_dt = 0.05
    interpolate = false
    max_residual_growth_factor = 1.01
  [../]
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
  picard_max_its = 2
  nl_abs_tol = 1e-09
  nl_rel_tol = 1e-06
  picard_abs_tol = 1e-08
  picard_rel_tol = 1e-04
[]

[Outputs]
  exodus = true
[]

[MultiApps]
  [./sub]
    type = TransientMultiApp
    app_type = MooseTestApp
    positions = '0 0 0'
    input_files = picard_sub.i
  [../]
[]

[Transfers]
  [./v_from_sub]
    type = MultiAppNearestNodeTransfer
    direction = from_multiapp
    multi_app = sub
    source_variable = v
    variable = v
  [../]
[]
