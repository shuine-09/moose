[Mesh]
  type = GeneratedMesh
  dim = 2
  xmin = -0.5
  xmax = 0.5
  ymin = -0.5
  ymax = 0.5
  nx = 32
  ny = 32
  elem_type = QUAD9
[]

[Variables]
  [./phi]
    order = SECOND
  [../]
[]

[AuxVariables]
  [./phi_0]
    order = SECOND
  [../]
  # [./marker]
  #   family = MONOMIAL
  #   order = CONSTANT
  # [../]
[]

[Kernels]
  [./time]
    type = TimeDerivative
    variable = phi
  [../]
  [./reinit]
    type = LevelSetOlssonReinitialization
    variable = phi
    phi_0 = phi_0
    epsilon = 0.04 #0.022097
  [../]
[]

[Problem]
  type = LevelSetReinitializationProblem
[]

[UserObjects]
  [./arnold]
    type = LevelSetOlssonTerminator
    tol = 0.5
    min_steps = 3
  [../]
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  start_time = 0
  num_steps = 5
  nl_abs_tol = 1e-14
  scheme = crank-nicolson
  line_search = none
  petsc_options_iname = '-pc_type -pc_sub_type'
  petsc_options_value = 'asm      ilu'
  dt = 0.01
[]

[Outputs]
[]
