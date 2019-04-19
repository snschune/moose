[Problem]
  kernel_coverage_check = false
[]

[Mesh]
  type = FileMesh
  file = generate_out.e
[]

[Variables]
  [./d]
  [../]
[]

[AuxVariables]
  [./restricted]
    family = MONOMIAL
    order = CONSTANT
    #block = 1
    block = 800
    initial_from_file_var = restricted
  [../]
[]

[Postprocessors]
  [./restricted]
    type = ElementAverageValue
    variable = restricted
    #block = 1
    block = 800
  [../]
[]

[Executioner]
  type = Steady
[]
