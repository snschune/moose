[Mesh]
  file = 2dcontact_collide.e
[]

[Variables]
  [./u]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxVariables]
  [./penetration]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  [./diff]
    type = Diffusion
    variable = u
  [../]
[]

[AuxBCs]
  [./penetration]
    type = PenetrationAux
    variable = penetration
    boundary = 2
    paired_boundary = 3
  [../]
[]

[BCs]
  [./block1_left]
    type = DirichletBC
    variable = u
    boundary = 1
    value = 0
  [../]
  [./block1_right]
    type = DirichletBC
    variable = u
    boundary = 2
    value = 1
  [../]
  [./block2_left]
    type = DirichletBC
    variable = u
    boundary = 3
    value = 0
  [../]
  [./block2_right]
    type = DirichletBC
    variable = u
    boundary = 4
    value = 1
  [../]
[]

[Executioner]
  type = Steady

  # Preconditioned JFNK (default)
  solve_type = 'PJFNK'

[]

[Outputs]
  output_initial = true
  exodus = true
  [./console]
    type = Console
    perf_log = true
  [../]
[]
