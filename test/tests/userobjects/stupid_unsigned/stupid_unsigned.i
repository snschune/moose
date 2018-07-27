[Mesh]
  type = GeneratedMesh
  dim = 2
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1
  nx = 5
  ny = 5
  elem_type = QUAD4
[]

[Problem]
  kernel_coverage_check = false
[]

[UserObjects]
  [./ud]
    type = StupidReadUnsignedInt
    dump_unsigned = -1
  [../]
[]

[Variables]
  [./u]
    family = LAGRANGE
    order = FIRST
  [../]
[]

[Executioner]
  type = Steady
[]
