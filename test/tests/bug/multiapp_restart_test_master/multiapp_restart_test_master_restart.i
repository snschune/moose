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
  type = Transient

  scheme = implicit-euler

  dt = 0.1
  end_time = 1.0

  l_tol = 1e-4
  nl_rel_tol = 1e-7
  nl_abs_tol = 1e-8
  l_max_its = 75
  nl_max_its = 50

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -ksp_gmres_restart'
  petsc_options_value = 'lu superlu_dist 75'

[]

[Outputs]
  exodus = true
[]

[MultiApps]
  [./sub]
    type = TransientMultiApp
    input_files = multiapp_restart_test_sub_restart.i
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
