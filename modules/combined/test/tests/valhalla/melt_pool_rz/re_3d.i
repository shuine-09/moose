[Mesh]
  type = GeneratedMesh
  dim = 3
  xmin = -1
  xmax = 1
  ymin = -1
  ymax = 1
  zmin = -1
  zmax = 1
  nx = 50
  ny = 50
  nz = 50
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
    epsilon = 0.2
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
