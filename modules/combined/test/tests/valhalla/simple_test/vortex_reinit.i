[Mesh]
  type = GeneratedMesh
  dim = 2
  xmax = 1
  ymax = 1
  nx = 16
  ny = 16
  uniform_refine = 2
  elem_type = QUAD9
  second_order = true
[]

[AuxVariables]
  [./velocity]
    family = LAGRANGE_VEC
  [../]
[]

[AuxKernels]
  [./vec]
    type = VectorFunctionAux
    variable = velocity
    function = velocity_func
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
[]

[Variables]
  [phi]
    family = LAGRANGE
  []
  [grad_ls]
    family = LAGRANGE_VEC
  []
[]

[Functions]
  [phi_exact]
    type = LevelSetOlssonBubble
    epsilon = 0.03
    center = '0.5 0.75 0'
    radius = 0.15
  []
  [./velocity_func]
    type = LevelSetOlssonVortex
    reverse_time = 2
  [../]
[]

[ICs]
  [phi_ic]
    type = FunctionIC
    function = phi_exact
    variable = phi
  []
[]

[Kernels]
  [time]
    type = TimeDerivative
    variable = phi
  []
  [advection]
    type = LevelSetAdvection
    velocity = velocity
    variable = phi
  []
  [advection_supg]
    type = LevelSetAdvectionSUPG
    velocity = velocity
    variable = phi
  []
  [time_supg]
    type = LevelSetTimeDerivativeSUPG
    velocity = velocity
    variable = phi
  []

  [grad_ls]
    type = VariableGradientRegularization
    regularized_var = phi
    variable = grad_ls
  []

  [reinit]
    type = LevelSetGradientRegularizationReinitialization
    variable = phi
    level_set_gradient = grad_ls
    epsilon = 0.03
    gamma = 1
  []
  [reinit_supg]
    type = LevelSetGradientRegularizationReinitializationSUPG
    variable = phi
    level_set_gradient = grad_ls
    epsilon = 0.03
    gamma = 1
    velocity = velocity
  []
[]

[Postprocessors]
  [area]
    type = LevelSetVolume
    threshold = 0.5
    variable = phi
    location = outside
    execute_on = 'initial timestep_end'
  []
  [cfl]
    type = LevelSetCFLCondition
    velocity = velocity
    execute_on = 'initial timestep_end'
  []
[]

[Problem]
  type = LevelSetProblem
[]

[Preconditioning/smp]
    type = SMP
    full = true
[]

[Executioner]
  type = Transient
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -ksp_type'
  petsc_options_value = 'lu superlu_dist preonly'
  line_search = 'none'
  solve_type = NEWTON
  start_time = 0
  end_time = 2
  #scheme = crank-nicolson
  [TimeStepper]
    type = PostprocessorDT
    postprocessor = cfl
    scale = 0.8
  []
[]

[Outputs]
  csv = true
  exodus = true
[]
