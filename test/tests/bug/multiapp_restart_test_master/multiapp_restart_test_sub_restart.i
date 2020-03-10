[Mesh]
  file = multiapp_restart_test_master_ss_out_sub0_Checkpoint_cp/LATEST
[]

[Problem]
  restart_file_base = multiapp_restart_test_master_ss_out_sub0_Checkpoint_cp/LATEST
  skip_additional_restart_data = true
  force_restart = true
[]

[Variables]
  [./u]
  [../]
[]

[Functions]
  [./conditional_function]
     type = ParsedFunction
     value = 't >= 0.5'
     direction = left
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
  [./trip_valve]
    type = FunctionValuePostprocessor
    function = conditional_function
    execute_on = 'INITIAL LINEAR TIMESTEP_END FINAL'
  [../]
  [./time_step_pp]
    type = TimestepSize
    # allow_duplicate_execution_on_initial = true
    execute_on = 'INITIAL LINEAR NONLINEAR TIMESTEP_END FINAL'
  [../]
[]

[Executioner]
  type = Transient # Pseudo transient to reach steady state.
  solve_type = 'PJFNK'
  petsc_options = ' -snes_converged_reason '
  line_search = l2

  # Problem time parameters.
  dt = 1e+15 # Let the master app control time steps.
  reset_dt = true
  start_time = 0.0

  # Iterations parameters.
  l_max_its = 50
  l_tol     = 1e-6

  nl_max_its = 25
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-9
[]

[Outputs]
  exodus = true
  csv = true
[]
