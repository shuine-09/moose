[Mesh]
  type = GeneratedMesh
  dim = 2
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 2
  nx = 32
  ny = 64
  elem_type = QUAD4
  uniform_refine = 1
[]

[Variables]
  [./phi]
    order = FIRST
  [../]
  [./normal]
    family = LAGRANGE_VEC
  [../]
[]

[AuxVariables]
  [./phi_0]
    order = FIRST
  [../]
  # [./marker]
  #   family = MONOMIAL
  #   order = CONSTANT
  # [../]
[]

[Kernels]
  [./grad_phi]
    type = LevelSetGradientRegulization
    c = phi_0
    variable = normal
    normal_regulization = true
  [../]
  [./time]
    type = TimeDerivative
    variable = phi
  [../]
  [./reinit]
    type = LevelSetOlssonReinitialization
    variable = phi
    phi_0 = phi_0
    epsilon = 0.02 #0.022097
    grad_c = normal
  [../]
[]

[Problem]
  type = LevelSetReinitializationProblem
[]

[UserObjects]
  [./arnold]
    type = LevelSetOlssonTerminator
    tol = 0.1
    min_steps = 10
  [../]
[]

[Postprocessors]
  [./area]
  type = LevelSetVolume
  threshold = 0.5
  variable = phi
  location = outside
  execute_on = 'timestep_end'
[../]
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  start_time = 0
  num_steps = 100
  nl_abs_tol = 1e-14
  scheme = crank-nicolson
  line_search = none
  petsc_options_iname = '-pc_type -pc_sub_type'
  petsc_options_value = 'asm      ilu'
  dt = 0.001
  # [./TimeStepper]
  #   type = IterationAdaptiveDT
  #   dt = 0.001
  #   optimal_iterations = 5
  #   growth_factor = 5
  # [../]
[]

[Outputs]
  exodus = true
  execute_on = 'TIMESTEP_END'
[]
