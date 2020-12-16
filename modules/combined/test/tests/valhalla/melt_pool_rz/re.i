[Mesh]
  type = GeneratedMesh
  dim = 2
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1
  nx = 20
  ny = 20
  uniform_refine = 3
[]

[Problem]
  coord_type = RZ
[]

[Variables]
  [./phi]
  [../]
[]

[AuxVariables]
  [./phi_0]
  [../]
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
    epsilon = 0.01
  [../]
[]

[Problem]
  type = LevelSetReinitializationProblem
[]

[UserObjects]
  [./arnold]
    type = LevelSetOlssonTerminator
    tol = 1
  [../]
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  start_time = 0
  num_steps = 100
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  scheme = crank-nicolson
  petsc_options_iname = '-pc_type -pc_sub_type -ksp_gmres_restart'
  petsc_options_value = 'hypre    boomeramg    300'
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.001
    optimal_iterations = 5
    growth_factor = 5
  [../]
[]

[Outputs]
[]
