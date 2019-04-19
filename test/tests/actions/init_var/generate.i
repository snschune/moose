[MeshGenerators]
  [./cartesian]
    type = CartesianMeshGenerator
    dim = 2
    dx = '1 1'
    dy = '1 1'
    subdomain_id = '1 1 1 800'
  [../]
[]

[Problem]
  kernel_coverage_check = false
[]

[Mesh]
  type = MeshGeneratorMesh
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
    initial_condition = 4
  [../]
[]

[Executioner]
  type = Steady
[]

[Postprocessors]
  [./restricted]
    type = ElementAverageValue
    variable = restricted
    #block = 1
    block = 800
  [../]
[]

[Outputs]
  exodus = true
[]
