[Mesh]
  type = FileMesh
  file = rectangle.e
  dim = 2
[]

[Variables]
  [./u]
  [../]
[]

[Functions]
  [./func]
    type = ParsedFunction
    value = x*y+t
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
[]

[BCs]
  [./left]
    type = DirichletBC
    variable = u
    boundary = 1
    value = 0
  [../]
  [./right]
    type = DirichletBC
    variable = u
    boundary = 2
    value = 1
  [../]
[]

[Materials]
  [./block]
    type = OutputTestMaterial
    block = '1 2'
    output_properties = tensor_property
    outputs = exodus
  [../]
  [./boundary_1]
    type = OutputTestMaterial
    boundary = 1
    output_properties = real_property
    outputs = exodus
    real_factor = 2
  [../]
  [./boundary_2]
    type = OutputTestMaterial
    boundary = 2
    output_properties = 'real_property vector_property'
    real_factor = 2
    outputs = exodus
  [../]
[]

[Executioner]
  # Preconditioned JFNK (default)
  type = Transient
  num_steps = 20
  dt = 0.1
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  output_initial = true
  exodus = true
  [./console]
    type = Console
    perf_log = true
    nonlinear_residuals = true
    linear_residuals = true
  [../]
[]