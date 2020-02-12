[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 0.01
    ymin = 0
    ymax = 0.01
    nx = 50
    ny = 50
    elem_type = QUAD4
  []
  [corner_node]
    type = ExtraNodesetGenerator
    new_boundary = 'pinned_node'
    coord = '0.0 0.01'
    input = gen
  []
[]

[Adaptivity]
  steps = 3
  marker = box
  max_h_level = 3
  initial_steps = 3
  stop_time = 1.0e-10
  [Markers]
    [box]
      bottom_left = '0.000 0.004 0'
      inside = refine
      top_right = '0.01 0.006 0'
      outside = do_nothing
      type = BoxMarker
    []
  []
[]

[ICs]
  [ls_ic]
    type = FunctionIC
    function = ls_exact
    variable = ls
  []
  [velocity]
    type = VectorConstantIC
    x_value = 1e-15
    y_value = 1e-15
    variable = velocity
  []
[]

[Variables]
  [ls]
    order = FIRST
  []
  [velocity]
    family = LAGRANGE_VEC
    order = FIRST
  []
  [p]
    order = FIRST
  []
  [curvature]
  []
  [temp]
    initial_condition = 300
  []
[]

[AuxVariables]
  [vel_x]
  []
  [vel_y]
  []
  [density]
    family = MONOMIAL
    order = CONSTANT
  []
  [specific_heat]
    family = MONOMIAL
    order = CONSTANT
  []
  [enthalpy]
    family = MONOMIAL
    order = CONSTANT
  []
  [thermal_conductivity]
    family = MONOMIAL
    order = CONSTANT
  []
  [liquid_mass_fraction]
    family = MONOMIAL
    order = CONSTANT
  []
  [solid_mass_fraction]
    family = MONOMIAL
    order = CONSTANT
  []
  [mu]
    family = MONOMIAL
    order = CONSTANT
  []
  [permeability]
    family = MONOMIAL
    order = CONSTANT
  []
  [melt_pool_momentum_source_x]
    family = MONOMIAL
    order = CONSTANT
  []
  [melt_pool_momentum_source_y]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[Functions]
  [ls_exact]
    type = LevelSetOlssonPlane
    epsilon = 0.00004
  []
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
  [density]
    type = MaterialRealAux
    variable = density
    property = rho
    execute_on = timestep_end
  []
  [specific_heat]
    type = MaterialRealAux
    variable = specific_heat
    property = specific_heat
    execute_on = timestep_end
  []
  [enthalpy]
    type = MaterialRealAux
    variable = enthalpy
    property = enthalpy
    execute_on = timestep_end
  []
  [thermal_conductivity]
    type = MaterialRealAux
    variable = thermal_conductivity
    property = thermal_conductivity
    execute_on = timestep_end
  []
  [liquid_mass_fraction]
    type = MaterialRealAux
    variable = liquid_mass_fraction
    property = liquid_mass_fraction
    execute_on = timestep_end
  []
  [solid_mass_fraction]
    type = MaterialRealAux
    variable = solid_mass_fraction
    property = solid_mass_fraction
    execute_on = timestep_end
  []
  [viscosity]
    type = MaterialRealAux
    variable = mu
    property = mu
    execute_on = timestep_end
  []
  [permeability]
    type = MaterialRealAux
    variable = permeability
    property = permeability
    execute_on = timestep_end
  []
  [melt_pool_momentum_source_x]
    type = MaterialRealVectorValueAux
    variable = melt_pool_momentum_source_x
    property = melt_pool_momentum_source
    component = 0
    execute_on = timestep_end
  []
  [melt_pool_momentum_source_y]
    type = MaterialRealVectorValueAux
    variable = melt_pool_momentum_source_y
    property = melt_pool_momentum_source
    component = 1
    execute_on = timestep_end
  []
[]

[Kernels]
  [curvature]
    type = LevelSetCurvatureRegulization
    level_set = ls
    variable = curvature
    epsilon = 2e-4
  []

  [time]
    type = TimeDerivative
    variable = ls
  []

  # [advection_supg]
  #   type = LevelSetAdvectionSUPG
  #   velocity_x = vel_x
  #   velocity_y = vel_y
  #   variable = ls
  # []
  #
  # [time_supg]
  #   type = LevelSetTimeDerivativeSUPG
  #   velocity_x = vel_x
  #   velocity_y = vel_y
  #   variable = ls
  # []
  #
  # [advection]
  #   type = LevelSetAdvection
  #   velocity_x = vel_x
  #   velocity_y = vel_y
  #   variable = ls
  # []
  #
  # [mass_addition]
  #   type = LevelSetPowderMass
  #   variable = ls
  #   mass_rate = 40.0e-4
  #   mass_radius = 0.25e-3
  #   laser_location = location
  # []

  [heat_time]
    type = ADHeatConductionTimeDerivative
    specific_heat = specific_heat
    density_name = rho
    variable = temp
  []

  [heat_conv]
    type = HeatConvection
    variable = temp
    velocity = velocity
    density = 'rho'
    specific_heat = 'specific_heat'
  []

  [heat_cond]
    type = ADHeatConduction
    thermal_conductivity = thermal_conductivity
    variable = temp
  []

  [heat_source]
    type = MeltPoolHeatSource
    variable = temp
    level_set = ls
    laser_power = 75
    effective_beam_radius = 0.25e-3
    absorption_coefficient = 0.27
    heat_transfer_coefficient = 100
    StefanBoltzmann_constant = 5.67e-8
    material_emissivity = 0.59
    ambient_temperature = 300
    laser_location = location
  []

  [mass]
    type = INSADMass
    variable = p
  []
  [mass_pspg]
    type = INSADMassPSPG
    variable = p
  []

  [momentum_time]
    type = INSADMomentumTimeDerivative
    variable = velocity
  []

  [momentum_convection]
    type = INSADMomentumAdvection
    variable = velocity
  []

  [momentum_viscous]
    type = INSADMomentumViscous
    variable = velocity
  []

  [momentum_pressure]
    type = INSADMomentumPressure
    variable = velocity
    p = p
    integrate_p_by_parts = true
  []

  [momentum_source]
    type = MeltPoolMomentumSource
    variable = velocity
  []

  [momentum_supg]
    type = INSADMomentumSUPG
    variable = velocity
    velocity = velocity
  []
[]

[UserObjects]
  [location]
    type = MeltPoolLevelSetLocation
    laser_center_x = '0.003-3.0e-04+6*1e-3*t'
    laser_center_y = '0.005'
    level_set = ls
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
[]

[Materials]
  [thermal]
    type = LevelSetThermalMaterial
    level_set = ls
    temperature = temp
    c_g = 300
    c_s = 500
    c_l = 500
    k_g = 0.017
    k_s = 31.8724
    k_l = 209.3
    solidus_temperature = 1648
    latent_heat = 2.5e5
  []
  [mushy]
    type = MushyZoneMaterial
    temperature = temp
    liquidus_temperature = 1673
    solidus_temperature = 1648
    rho_s = 8000
    rho_l = 8000
  []
  [melt]
    type = LevelSetMeltMaterial
    level_set = ls
    rho_g = 1.184
    rho_s = 8000
    rho_l = 8000
    mu_g = 1.81e-5
    mu_l = 0.1
    mu_s = 100
    permeability_constant = 1e-8
  []
  [ins_mat]
    type = MeltPoolLevelSetINSMaterial
    velocity = velocity
    pressure = p
    transient_term = true
    integrate_p_by_parts = true
    alpha = .1
    level_set = ls
    temperature = temp
    curvature = curvature
    surface_tension = 1.169
    thermal_capillary = -4.3e-4
    thermal_expansion = 1.45e-6
    reference_temperature = 300
    rho_l = 8000
  []
[]

[BCs]
  [no_slip]
    type = ADVectorFunctionDirichletBC
    variable = velocity
    boundary = '1 3 2 4'
  []
  [pressure_pin]
    type = DirichletBC
    variable = p
    boundary = 'pinned_node'
    value = 0
  []
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
    solve_type = 'NEWTON'
  []
[]

# [MultiApps]
#   [reinit]
#     type = LevelSetReinitializationMultiApp
#     input_files = 'ls_reinit.i'
#     execute_on = TIMESTEP_END
#   []
# []
#
#
#
# [Transfers]
#   # [marker_to_sub]
#   #   type = LevelSetMeshRefinementTransfer
#   #   multi_app = reinit
#   #   source_variable = marker
#   #   variable = marker
#   #   check_multiapp_execute_on = false
#   # []
#   [to_sub]
#     type = MultiAppCopyTransfer
#     source_variable = ls
#     variable = ls
#     direction = to_multiapp
#     multi_app = reinit
#     execute_on = 'timestep_end'
#   []
#
#   [to_sub_init]
#     type = MultiAppCopyTransfer
#     source_variable = ls
#     variable = ls_0
#     direction = to_multiapp
#     multi_app = reinit
#     execute_on = 'timestep_end'
#   []
#
#   [from_sub]
#     type = MultiAppCopyTransfer
#     source_variable = ls
#     variable = ls
#     direction = from_multiapp
#     multi_app = reinit
#     execute_on = 'timestep_end'
#   []
# []

[Executioner]
  type = Transient
  solve_type = NEWTON
  start_time = 0
  dt = 1e-3
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type -ksp_type'
  petsc_options_value = 'lu superlu_dist preonly'
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6
  nl_max_its = 15
  l_tol = 1e-6
  l_max_its = 20
  #dtmin = 1e-3
  dtmax = 1
  end_time = 1
  line_search = 'none'
  automatic_scaling = true
  skip_exception_check = true
  nl_div_tol = 1e10
[]

[Outputs]
  exodus = true
[]
