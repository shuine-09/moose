[Mesh]
  # type = GeneratedMesh
  # dim = 2
  # xmax = 1
  # ymax = 1
  # nx = 16
  # ny = 16
  # uniform_refine = 3
  # elem_type = QUAD4
  # second_order = true
  type = FileMesh
  file = 'tri_mesh.e'
[]

[Variables/phi]
    family = LAGRANGE
[]

[AuxVariables]
  [phi_0]
    family = LAGRANGE
  []
  [marker]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[Kernels]
  [time]
    type = TimeDerivative
    variable = phi
  []
  [reinit]
    type = LevelSetOlssonReinitialization
    variable = phi
    phi_0 = phi_0
    epsilon = 0.0063
  []
[]

[Problem]
  type = LevelSetReinitializationProblem
[]

[UserObjects]
  [arnold]
    type = LevelSetOlssonTerminator
    tol = 0.1
    min_steps = 10
  []
[]

[Preconditioning/smp]
    type = SMP
    full = true
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  start_time = 0
  num_steps = 100
  nl_abs_tol = 1e-14
  scheme = crank-nicolson
  line_search = none
  dt = 0.01
[]

[Outputs]
[]
