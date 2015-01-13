[Mesh]
 type =  GeneratedMesh
 dim = 2
 xmin = 0
 xmax = 4
 ymin = 0
 ymax = 4
 elem_type = QUAD4
 nx = 8
 ny = 8
 uniform_refine = 0
[]

[Variables]
  active = 'u'

  [./u]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  active = 'diff'

  [./diff]
    type = Diffusion
    variable = u
  [../]
[]

[BCs]
  active = 'left right'

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
  [./sigma_total]
    type = GenericConstantMaterial
    block = 0
    prop_names = 'sigma_t'
    prop_values = '1.0'
  [../]
  [./mult_mat1]
    type = GenericConstantMaterial
    block = 0
    prop_names = 'mult1 mult2 mult3'
    prop_values = '2.3 -1.0 0.32'
  [../]
[]

[Executioner]
  type = Steady

  solve_type = 'NEWTON'
[]

[Outputs]
  file_base = out
  exodus = true
  output_on = 'initial timestep_end'
  [./console]
    type = Console
    perf_log = true
    output_on = 'timestep_end failed nonlinear linear'
  [../]
[]
