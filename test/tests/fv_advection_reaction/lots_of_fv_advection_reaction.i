[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 10
[]

[Problem]
  kernel_coverage_check = false
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  l_max_its = 100
  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = ' lu       NONZERO               1e-10'
[]

[Outputs]
  execute_on = 'timestep_end'
  exodus = true
[]

[Testing]
  [LotsOfFVAdvectionReaction]
    vel = '1 0 0'
    inlet_boundaries = left
    outlet_boundaries = right
    number = 10
  []
[]
