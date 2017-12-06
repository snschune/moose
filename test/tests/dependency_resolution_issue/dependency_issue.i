[Mesh]
  type = GeneratedMesh
  dim = 2
  xmin = 0
  xmax = 1
  nx = 2
  ymin = 0
  ymax = 1
  ny = 2
[]

[Variables]
  [./u]
    initial_condition = 1
  [../]
[]

[Kernels]
  [./time]
    type = TimeDerivative
    variable = u
  [../]

  [./diff]
    type = Diffusion
    variable = u
  [../]
[]

[AuxVariables]
  [./scaled_u]
  [../]

  [./dummy]
  [../]
[]

[AuxKernels]
  [./normal_aux]
    type = NormalizationAux
    source_variable = u
    variable = scaled_u
    normalization = scale1
    execute_on = timestep_end
  [../]

  [./dummy]
    type = NormalizationAux
    source_variable = u
    variable = dummy
    normalization = scale2
    execute_on = timestep_begin
  [../]
[]

[Postprocessors]
  [./total_u]
    type = ElementIntegralVariablePostprocessor
    variable = u
    execute_on = timestep_begin
  [../]

  # scale 2 is always in preaux and it should be because of the dummy AuxKernel
  # correctly putting it there.
  # scale 1 should be in postaux because no AuxKernel in timestep_begin depends
  # on it, but it is bumped to preaux by normal_aux [note: execute_on = timestep_end]
  # you can check that easily by looking at output:
  #
  #  scale 1 value at time = 1 should be equal to total_u because GeneralUO executes
  #  after ElementUO. This shows scale1 was bumped to pre_aux.
  #  +----------------+----------------+----------------+----------------+
  #  | time           | scale1         | scale2         | total_u        |
  #  +----------------+----------------+----------------+----------------+
  #  |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
  #  |   1.000000e+00 |   0.000000e+00 |   0.000000e+00 |   1.000000e+00 |
  #  |   2.000000e+00 |   1.000000e+00 |   1.000000e+00 |   2.260024e-01 |
  #  +----------------+----------------+----------------+----------------+
  #
  #  we can check that very easily by changing AuxKernels/normal_aux/normalization=scale2
  #  giving the following answer
  #  +----------------+----------------+----------------+----------------+
  #  | time           | scale1         | scale2         | total_u        |
  #  +----------------+----------------+----------------+----------------+
  #  |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |   0.000000e+00 |
  #  |   1.000000e+00 |   1.000000e+00 |   0.000000e+00 |   1.000000e+00 |
  #  |   2.000000e+00 |   2.260024e-01 |   1.000000e+00 |   2.260024e-01 |
  #  +----------------+----------------+----------------+----------------+
  #
  [./scale1]
    type = ScalePostprocessor
    value = total_u
    scaling_factor = 1
    execute_on = timestep_begin
  [../]

  [./scale2]
    type = ScalePostprocessor
    value = total_u
    scaling_factor = 1
    execute_on = timestep_begin
  [../]
[]

[BCs]
  [./left]
    type = DirichletBC
    variable = u
    boundary = 1
    value = 0
  [../]
[]

[Executioner]
  type = Transient
  dt = 1.0
  end_time = 2.0
[]
