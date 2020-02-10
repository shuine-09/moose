[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = -0.5
    xmax = 0.5
    ymin = -0.5
    ymax = 0.5
    nx = 32
    ny = 32
    elem_type = QUAD9
  []
  [./corner_node]
    type = ExtraNodesetGenerator
    new_boundary = 'pinned_node'
    nodes = '0'
    input = gen
  [../]
[]

# [GlobalParams]
#   gravity = '0 -0.1 0'
# []

[AuxVariables]
  [vel_x]
  []
  [vel_y]
  []
[]

[ICs]
  [velocity]
    type = VectorConstantIC
    x_value = 1e-15
    y_value = 1e-15
    variable = velocity
  []
  [./phi_ic]
    type = FunctionIC
    function = phi_exact
    variable = phi
  [../]
[]

[AuxKernels]
  [vel_x]
    type = VectorVariableComponentAux
    variable = vel_x
    vector_variable = velocity
    component = 'x'
  []
  [vel_y]
    type = VectorVariableComponentAux
    variable = vel_y
    vector_variable = velocity
    component = 'y'
  []
[]

[Variables]
  [./velocity]
    family = LAGRANGE_VEC
  [../]
  [./grad_phi]
    family = LAGRANGE_VEC
  [../]
  [./curvature]
  [../]
  [./p]
  [../]
  [./phi]
    order = SECOND
  [../]
[]

[Problem]
  type = LevelSetProblem
[]

[Functions]
  [./phi_exact]
    type = LevelSetOlssonBubble
    epsilon = 0.04
    center = '0.5 0.5 0'
    radius = 0.25
  [../]
[]

[Kernels]
  [./grad_phi]
    type = LevelSetGradientRegulization
    c = phi
    variable = grad_phi
  [../]
  [./curvature]
    type = LevelSetCurvatureRegulization
    grad_c = grad_phi
    variable = curvature
    epsilon = 0.005
  [../]
  [./momentum_time]
    type = INSADMomentumTimeDerivative
    variable = velocity
  [../]
  [./mass]
    type = INSADMass
    variable = p
  [../]
  [./mass_pspg]
    type = INSADMassPSPG
    variable = p
  [../]

  [./momentum_convection]
    type = INSADMomentumAdvection
    variable = velocity
  [../]

  [./momentum_viscous]
    type = INSADMomentumViscous
    variable = velocity
  [../]

  [./momentum_pressure]
    type = INSADMomentumPressure
    variable = velocity
    p = p
    integrate_p_by_parts = true
  [../]

  [./momentum_supg]
    type = INSADMomentumSUPG
    variable = velocity
    velocity = velocity
  [../]

  [./surface_tension]
    type = INSTwoPhaseSurfaceTension
    variable = velocity
    grad_c = grad_phi
    curvature = curvature
    sigma = 0.1
  [../]

  [./time]
    type = TimeDerivative
    variable = phi
  [../]
  [./advection]
    type = LevelSetAdvection
    velocity_x = vel_x
    velocity_y = vel_y
    variable = phi
  [../]
  [./advection_supg]
    type = LevelSetAdvectionSUPG
    velocity_x = vel_x
    velocity_y = vel_y
    variable = phi
  [../]
  [./time_supg]
    type = LevelSetTimeDerivativeSUPG
    velocity_x = vel_x
    velocity_y = vel_y
    variable = phi
  [../]
[]

[BCs]
  [./no_slip]
    type = VectorFunctionDirichletBC
    variable = velocity
    boundary = 'bottom right left top'
  [../]

  [./pressure_pin]
    type = DirichletBC
    variable = p
    boundary = 'pinned_node'
    value = 0
  [../]
[]

[Materials]
  [rho]
    type = LevelSetTwoPhaseMaterial
    prop_name = rho
    prop_value_1 = 10000
    prop_value_2 = 10000
    phi = phi
    outputs = exodus
  []
  [mu]
    type = LevelSetTwoPhaseMaterial
    prop_name = mu
    prop_value_1 = 1
    prop_value_2 = 1
    phi = phi
    outputs = exodus
  []
  # [./const]
  #   type = GenericConstantMaterial
  #   prop_names = 'rho mu'
  #   prop_values = '1  1'
  # [../]
  [ins_mat]
    type = INSADTauMaterial
    velocity = velocity
    pressure = p
    transient_term = true
    integrate_p_by_parts = true
    alpha = .1
  []
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
    solve_type = 'NEWTON'
  [../]
[]

[Postprocessors]
  [./cfl]
    type = LevelSetCFLCondition
    velocity_x = vel_x
    velocity_y = vel_y
    execute_on = 'initial timestep_end'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  start_time = 0
  dt = 1
  end_time = 1000
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu mumps'
  line_search = 'none'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-10
  nl_max_its = 15
  l_tol = 1e-6
  l_max_its = 20
  dtmax = 1
  [./TimeStepper]
    type = PostprocessorDT
    postprocessor = cfl
    scale = 0.5
    dt = 1
  [../]
[]

[Outputs]
  exodus = true
[]

# [MultiApps]
#   [./reinit]
#     type = LevelSetReinitializationMultiApp
#     input_files = 'static_bubble_reinit_sub.i'
#     execute_on = TIMESTEP_END
#   [../]
# []
# 
# [Transfers]
#   # [./marker_to_sub]
#   #   type = LevelSetMeshRefinementTransfer
#   #   multi_app = reinit
#   #   source_variable = marker
#   #   variable = marker
#   #   check_multiapp_execute_on = false
#   # [../]
#   [./to_sub]
#     type = MultiAppCopyTransfer
#     source_variable = phi
#     variable = phi
#     direction = to_multiapp
#     multi_app = reinit
#     execute_on = 'timestep_end'
#   [../]
#
#   [./to_sub_init]
#     type = MultiAppCopyTransfer
#     source_variable = phi
#     variable = phi_0
#     direction = to_multiapp
#     multi_app = reinit
#     execute_on = 'timestep_end'
#   [../]
#
#   [./from_sub]
#     type = MultiAppCopyTransfer
#     source_variable = phi
#     variable = phi
#     direction = from_multiapp
#     multi_app = reinit
#     execute_on = 'timestep_end'
#   [../]
# []
