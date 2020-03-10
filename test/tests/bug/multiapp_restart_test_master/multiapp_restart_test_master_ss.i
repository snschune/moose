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
  [./from_subapp]
    initial_condition = -1000
  [../]
[]

[Kernels]
  [./diff]
    type = Diffusion
    variable = u
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

[Executioner]
  type = Steady
  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'

  picard_abs_tol = 1e-8
  picard_rel_tol = 1e-7
  picard_max_its = 2
  disable_picard_residual_norm_check = true
[]

[Outputs]
  exodus = true
[]

[MultiApps]
  [./sub]
    type = FullSolveMultiApp
    input_files = multiapp_restart_test_sub_ss.i
    keep_solution_during_restore = true
    max_procs_per_app = 1
  [../]
[]

[Postprocessors]
  [./from_subapp_average]
    type = ElementAverageValue
    variable = from_subapp
  [../]
[]

[Transfers]
  [./sub_average]
    type = MultiAppPostprocessorInterpolationTransfer
    direction = from_multiapp
    multi_app = sub
    variable = from_subapp
    postprocessor = sub_average
  [../]
[]
