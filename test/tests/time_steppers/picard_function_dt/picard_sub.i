[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 2
  ny = 2
[]

[GlobalParams]
  implicit = false
[]

[Variables]
  [./v]
  [../]
[]

[Kernels]
  [./diff_v]
    type = Diffusion
    variable = v
  [../]
  [./time]
    type = TimeDerivative
    variable = v
  [../]
  [./forcing]
    type = BodyForce
    variable = v
  [../]
[]

[BCs]
  [./bc]
    type = DirichletBC
    variable = v
    boundary = 'left right'
    value = 0
  [../]
[]

[Executioner]
  # Preconditioned JFNK (default)
  type = Transient
  scheme = 'explicit-euler'
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
  nl_abs_tol = 1e-06
  nl_rel_tol = 1e-05
  nl_max_its = 100
[]

[Postprocessors]
  [./norm]
    type = ElementIntegralVariablePostprocessor
    variable = v
  [../]
[]

[Outputs]
  exodus = true
[]
