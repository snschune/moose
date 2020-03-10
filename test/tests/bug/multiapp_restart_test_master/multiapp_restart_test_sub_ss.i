[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 2
  ny = 2
[]

[Variables]
  [./u]
  [../]
[]

[AuxVariables]
[]

[Kernels]
  [./diff]
    type = Diffusion
    variable = u
  [../]
  [./force]
    type = CoupledForce
    variable = u
    v = 1.0
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
  [./sub_average]
    type = ElementAverageValue
    variable = u
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'

  # Problem time parameters.
  reset_dt = true
  start_time = 0.0

  end_time = 1e+6
  dtmin    = 1e-2
  dtmax    = 5e+4

  # Linear/nonlinear iterations.
  l_max_its = 50
  l_tol     = 1e-6

  nl_max_its = 25
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-9

  # Steady state detection.
  steady_state_detection = true
  steady_state_tolerance = 1e-10
  steady_state_start_time = 1.0

  # Time step control.
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt                 = 1e-2
    cutback_factor     = 0.5
    growth_factor      = 2.00
    optimal_iterations = 25
  [../]
[]

[Outputs]
  exodus = true
  [./Checkpoint]
    type = Checkpoint
    additional_execute_on = 'FINAL'
  [../]
[]
