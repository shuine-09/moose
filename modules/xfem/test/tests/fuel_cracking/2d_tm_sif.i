[GlobalParams]
  density = 10431.0 #95% TD (TD = 10980)
  order = FIRST
  family = LAGRANGE
  energy_per_fission = 3.2e-11  # J/fission (205 Mev)
  displacements = 'disp_x disp_y'
  volumetric_locking_correction = true
[]

# [Problem]
#   type = ReferenceResidualProblem
# #  solution_variables = 'disp_x disp_y strain_zz temp'
#   solution_variables = 'disp_x disp_y clad_strain_zz clad_strain_zz temp'
# #  reference_residual_variables = 'saved_x saved_y saved_z saved_t'
#   reference_residual_variables = 'saved_x saved_y saved_z saved_z saved_t'
# #  acceptable_iterations = 10
# #  acceptable_multiplier = 10
# []

[Mesh]
  file = 2d_fuel_clad_refine_1.e
  patch_size = 3
[]

[XFEM]
  qrule = volfrac
  output_cut_plane = true
[]

[ICs]
  [./threshold]
    type = RandomIC
    block = pellet_type_1
    min = 126.75e6
    max = 133.25e6
    variable = threshold
  [../]
[]

[UserObjects]
  [./xfem_marker_uo]
    type = XFEMMaterialTensorMarkerUserObject
    execute_on = timestep_end
    tensor = stress
    #quantity = MaxPrincipal
    quantity = hoop
    point1 = '0, 0, 0'
    point2 = '0, 0, 1'
    #average = false
    average = true
    threshold = 1.5e8 #1.3e8
    initiate_on_boundary = 8
    random_range = 0.1
    #weibull = xfem_weibull
    use_weibull = true
  [../]
  [./xfem_energy_release_rate]
    type = XFEMEnergyReleaseRate
    execute_on = timestep_end
    block = pellet_type_1
    use_weibull = true
    average_h = h
    radius_inner = 2e-4
    radius_outer = 4e-4
    critical_energy_release_rate = 200
  [../]
  [./xfem_max_hoop_stress]
    type = XFEMMaxHoopStress
    execute_on = timestep_end
    block = pellet_type_1
    disp_x = disp_x
    disp_y = disp_y
    temp = temp
    youngs_modulus = 2.e11
    poissons_ratio = .345
    thermal_expansion = 10e-6
    radius_inner = 2e-4
    radius_outer = 4e-4
    average_h = h
    intersecting_boundary = 8
 [../]
 [./pair_qps]
   type = XFEMElementPairQPProvider
 [../]
 [./manager]
   type = XFEMElemPairMaterialManager
   material_names = 'material1'
   element_pair_qps = pair_qps
 [../]
[]

[Variables]
  [./temp]
    initial_condition = 298.0
  [../]
  [./fuel_strain_zz]
    order = FIRST
    family = SCALAR
  [../]
  [./clad_strain_zz]
    order = FIRST
    family = SCALAR
  [../]
[]

[AuxVariables]
  [./threshold]
  [../]
  [./fission_rate]
    block = pellet_type_1
  [../]
  [./burnup]
    block = pellet_type_1
  [../]
  [./fast_neutron_flux]
  [../]
  [./fast_neutron_fluence]
  [../]
  [./hoop_stress]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./hydrostatic_stress]
    order = CONSTANT
    family = MONOMIAL
  [../]

  [./gap_conductance]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./pid]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./saved_x]
  [../]
  [./saved_y]
  [../]
  [./saved_z]
  [../]
  [./saved_t]
  [../]
  [./fuel_radius]
  [../]
  [./weibull]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Functions]
  [./power_history]
    type = PiecewiseLinear
    # scale_factor = 1.8936e+4 #Convert from W/m to W/m3 (1/pi/(.0041m)^2)
    x = '0 200 1.2e4 2.2e4 3.2e4 4.2e4'
    y = '0 0 2.5e4 0 2.5e4 0'
  [../]
  [./axial_peaking_factors]
    type = ParsedFunction
    value = 1
  [../]
  [./coolant_pressure_ramp]
    type = PiecewiseLinear
    x = '0 200'
    y = '6.537e-3 1'
    scale_factor = 15.5e6
  [../]
  [./q]                                 # this is for fuel_relocation
    type = CompositeFunction
    functions = 'power_history axial_peaking_factors'
  [../]

  [./clad_axial_pressure]
    type = ParsedFunction
    value = '(coolant_pressure * 7.0885128e-5 - plenum_pressure * 5.489116e-5) / 1.5991e-6'
    vars = 'plenum_pressure coolant_pressure'
    vals = 'plenum_pressure coolant_pressure'
  [../]
  [./fuel_axial_pressure]
    type = ParsedFunction
    value = plenum_pressure
    vars = plenum_pressure
    vals = plenum_pressure
  [../]
[]

[Constraints]
  [./dispx_constraint]
    type = XFEMResidualOpening
    use_displaced_mesh = false
    component = 0
    variable = disp_x
    disp_x = disp_x
    disp_y = disp_y
    residual_opening = 0.001
    alpha = 1.0e14
    manager = manager
  [../]
  [./dispy_constraint]
    type = XFEMResidualOpening
    use_displaced_mesh = false
    component = 1
    variable = disp_y
    disp_x = disp_x
    disp_y = disp_y
    residual_opening = 0.001
    alpha = 1.0e14
    manager = manager
  [../]
[]

[Modules]
  [./TensorMechanics]
    [./Master]
      [./gps_fuel]
        add_variables = true
        scalar_out_of_plane_strain = fuel_strain_zz
        out_of_plane_pressure = fuel_axial_pressure
        block = pellet_type_1
        planar_formulation = GENERALIZED_PLANE_STRAIN
        save_in = 'saved_x saved_y'
        strain = FINITE
        eigenstrain_names = 'fuel_thermal_strain fuel_volumetric_strain'
        generate_output = 'stress_xx stress_yy stress_zz vonmises_stress hydrostatic_stress firstinv_strain strain_zz'
      [../]
      [./gps_clad]
        add_variables = true
        scalar_out_of_plane_strain = clad_strain_zz
        out_of_plane_pressure = clad_axial_pressure
        planar_formulation = GENERALIZED_PLANE_STRAIN
        block = clad
        save_in = 'saved_x saved_y'
        strain = FINITE
        eigenstrain_names = 'clad_thermal_strain'
        generate_output = 'stress_xx stress_yy stress_zz vonmises_stress hydrostatic_stress firstinv_strain strain_zz'
      [../]
    [../]
  [../]
[]

[Kernels]
  [./heat]
    type = HeatConduction
    variable = temp
    save_in = saved_t
  [../]
  [./heat_ie]
    type = HeatConductionTimeDerivative
    variable = temp
    save_in = saved_t
  [../]
  [./heat_source_fuel]
    type = NeutronHeatSource
    variable = temp
    block = pellet_type_1
    fission_rate = fission_rate
    save_in = saved_t
  [../]
[]

[AuxKernels]
  [./fissionrate]
    type = FissionRateAux
    variable= fission_rate
    block = pellet_type_1
    value = 5.91742e+14
    function = q
  [../]

  [./burnup]
    type = BurnupAux
    variable = burnup
    block = pellet_type_1
    fission_rate = fission_rate
    molecular_weight = 0.270
  [../]

  [./fast_neutron_flux]
    type = FastNeutronFluxAux
    variable = fast_neutron_flux
    block = clad
    rod_ave_lin_pow = power_history
    axial_power_profile = axial_peaking_factors
    factor = 2.34e+13
    execute_on = timestep_begin
  [../]

  [./fast_neutron_fluence]
    type = FastNeutronFluenceAux
    variable = fast_neutron_fluence
    block = clad
    fast_neutron_flux = fast_neutron_flux
    execute_on = timestep_begin
  [../]

  [./hoop_stress]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = hoop_stress
    scalar_type = HoopStress
    point1 = '0 0 0'
    point2 = '0 0 1'
    execute_on = timestep_end
  [../]

  [./pid]
    type = ProcessorIDAux
    variable = pid
  [../]
  [./gap_cond]
    type = MaterialRealAux
    property = gap_conductance
    variable = gap_conductance
    boundary = 10
  [../]

  [./fuel_radius]
    type = Radius
    variable = fuel_radius
    boundary = 8
    point1 = '0 0 0'
    point2 = '0 0 1'
    execute_on = linear
  [../]
  [./weibull_output]
    type = MaterialRealAux
    variable = weibull
    property = weibull
  [../]
[]

[Contact]
  [./pellet_clad_mechanical]
    master = 5
    slave = 10
    system = Constraint
    penalty = 1e+12 #1e7
    model = frictionless
    tangential_tolerance = 5e-4
    normal_smoothing_distance = 0.1
    normalize_penalty = true
#    master_slave_jacobian = true
#    connected_slave_nodes_jacobian = true
#    non_displacement_variables_jacobian = false
  [../]
[]

[ThermalContact]
  [./pellet_clad_thermal]
    type = GapHeatTransferLWR
    variable = temp
    master = 5
    slave = 10
    jump_distance_model = KENNARD
    interaction_layer = 0
    plenum_pressure = plenum_pressure
    contact_pressure = contact_pressure
    roughness_coef = 3.2
    roughness_clad = 1e-6
    roughness_fuel = 2e-6
#    quadrature = true
    tangential_tolerance = 5e-4
    normal_smoothing_distance = 0.1
  [../]
[]

[BCs]

  [./no_x_clad]
    type = PresetBC
    variable = disp_x
    boundary = '10002'
    value = 0.0
  [../]

  [./no_y_clad]
    type = PresetBC
    variable = disp_y
    boundary = '10003'
    value = 0.0
  [../]

  [./no_x_fuel]
    type = PresetBC
    variable = disp_x
    boundary = '10000'
    value = 0.0
  [../]

  [./no_y_fuel]
    type = PresetBC
    variable = disp_y
    boundary = '10000 10001'
    value = 0.0
  [../]

  [./Pressure]
    # apply coolant pressure on clad outer walls
    [./coolantPressure]
      boundary = '2'
      function = coolant_pressure_ramp
    [../]
  [../]

  [./PlenumPressure]
    # apply plenum pressure on clad inner walls and pellet surfaces
    [./plenumPressure]
      boundary = 9
      #BWR change: fill pressure differs
      initial_pressure = 2e6
      startup_time = 0.0
      output_initial_moles = initial_moles
      R = 8.3143
      temperature = interior_temp
      volume = gas_volume
      output = plenum_pressure
    [../]
  [../]

  # convective boundary on clad outer surface
  [./convective_clad_surface]
    type = ConvectiveFluxBC
    boundary = '2'
    variable = temp
    #BWR change convection coefficient and temperature
    rate = 100e3         # convection coefficient (h, W/m2-K)
    initial = 298.0        # initial coolant temp
    final = 561.0          # saturated at 1050 psi
    duration = 200     # duration of initial power ramp
  [../]
[]

[Materials]
  [./fuel_thermal]
    type = ThermalFuel
    block = pellet_type_1
    temp = temp
    burnup = burnup
    model = 4
  [../]

  [./elastic_tensor]
    type = ComputeIsotropicElasticityTensor
    block = pellet_type_1
    poissons_ratio = 0.345
    youngs_modulus = 2.0e11
  [../]

  [./fuel_thermal_expansion]        # thermal expansion strain for UO2
    type = ComputeThermalExpansionEigenstrain
    block = pellet_type_1
    thermal_expansion_coeff = 10.0e-6
    temperature = temp
    stress_free_temperature = 298.0
    eigenstrain_name = fuel_thermal_strain
  [../]

  [./fuel_volumetric_swelling]      # free expansion strains (swelling and densification) for UO2 (bison kernel)
    type = UO2VolumetricSwellingEigenstrain
    gas_swelling_model_type = MATPRO
    block = pellet_type_1
    temperature = temp
    burnup = burnup
    initial_fuel_density = 10431.0
    eigenstrain_name = fuel_volumetric_strain
  [../]

  [./fuel_stress]
    type = ComputeFiniteStrainElasticStress
    block = pellet_type_1
  [../]

  [./clad_thermal]
    type = HeatConductionMaterial
    block = clad
    thermal_conductivity = 16.0
    specific_heat = 330.0
  [../]

  [./clad_elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    block = clad
    youngs_modulus = 7.5e10
    poissons_ratio = 0.3
  [../]

  [./clad_thermal_expansion]        # thermal expansion strain for UO2
    type = ComputeThermalExpansionEigenstrain
    block = clad
    thermal_expansion_coeff = 5.0e-6
    temperature = temp
    stress_free_temperature = 298.0
    eigenstrain_name = clad_thermal_strain
  [../]

  [./clad_stress]
    type = ComputeFiniteStrainElasticStress
    block = clad
  [../]

  [./clad_density]
    type = Density
    block = clad
    density = 6551.0
  [../]

  [./fuel_density]
    type = Density
    block = pellet_type_1
  [../]

  [./weilbull]
    type = XFEMWeibullMaterial
    weibull_modulus = 15
    specimen_volume = 5.281e-5 #5.281e-5 # 7.6977e-05
    specimen_material_property = 1.0
    block = pellet_type_1
  [../]

  [./material1]
    type = MaximumNormalSeparation
    compute = false
    disp_x = disp_x
    disp_y = disp_y
  [../]

[]

[Dampers]
  [./limitT]
    type = MaxIncrement
    variable = temp
    max_increment = 25.0
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
#    full = true
    off_diag_row    = 'disp_x disp_y'
    off_diag_column = 'disp_y disp_x'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'

  petsc_options = '-snes_ksp_ew'
#  petsc_options = '-snes_linesearch_monitor'
#  petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type -pc_hypre_boomeramg_max_iter -pc_hypre_strong_threshold'
#  petsc_options_value = '201                hypre    boomeramg      2                            0.7'

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -mat_superlu_dist_fact'
  petsc_options_value = 'lu       superlu_dist    SamePattern'

  line_search = 'none'

  l_max_its = 25
  nl_max_its = 50
  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-8
#  l_tol = 1e-2

  start_time = 0.0
  dt = 100
  end_time = 40200

  # [./TimeStepper]
  #   type = IterationAdaptiveDT
  #   dt = 1e2
  #   optimal_iterations = 30
  #   iteration_window = 4
  #   timestep_limiting_function = power_history
  #   force_step_every_function_point = true
  # [../]


  [./Quadrature]
     order = third #fifth
     side_order = fifth #seventh
  [../]

  max_xfem_update = 10000
  verbose = true
[]

[Postprocessors]
  [./react_z]
    type = MaterialTensorIntegral
    rank_two_tensor = stress
    index_i = 2
    index_j = 2
  [../]
  [./clad_inner_vol]
   type = InternalVolume
    boundary = 7
  [../]

  [./pellet_volume]
    type = InternalVolume
    boundary = 8
  [../]

  [./gas_volume]
    type = InternalVolume
    boundary = 9
    execute_on = 'initial timestep_end'
    addition = 6.73395e-7 #pi*D^2/4 *h = pi*(8.36e-3)^2*(38.5405e-3)/4
  [../]

  [./interior_temp]
    type = SideAverageValue
    boundary = 9 #kg was 7
    variable = temp
    execute_on = 'initial timestep_end'
  [../]

  [./coolant_pressure]
    type = FunctionValuePostprocessor
    function = coolant_pressure_ramp
    execute_on = 'initial linear'
  [../]

  [./power_history]
    type = FunctionValuePostprocessor
    function = power_history
  [../]

  [./flux_from_clad]
    type = SideFluxIntegral
    variable = temp
    boundary = 5
    diffusivity = thermal_conductivity
  [../]

  [./flux_from_fuel]
    type = SideFluxIntegral
    variable = temp
    boundary = 10
    diffusivity = thermal_conductivity
  [../]

  [./dt]
    type = TimestepSize
  [../]

  [./residual]
    type = Residual
  [../]

  [./nl_its]
    type = NumNonlinearIterations
  [../]

  [./lin_its]
    type = NumLinearIterations
  [../]

  [./average_burnup]
    type = ElementAverageValue
    block = pellet_type_1
    variable = burnup
  [../]
  [./average_fissionrate]
    type = ElementAverageValue
    block = pellet_type_1
    variable = fission_rate
  [../]

  [./rod_total_power]
    type = ElementIntegralPower
    variable = temp
    fission_rate = fission_rate
    block = pellet_type_1
  [../]

  [./rod_input_power]
    type = FunctionValuePostprocessor
    function = power_history
    scale_factor = 0.05175  #BWR change: length of fuel stack in meters (5*pellet height)
  [../]

  [./min_strain_zz]
    type = ElementExtremeValue
    variable = strain_zz
    value_type = min
  [../]

  [./max_strain_zz]
    type = ElementExtremeValue
    variable = strain_zz
    value_type = max
  [../]

  [./average_fuel_radius]
    type = AverageNodalVariableValue
    variable = fuel_radius
  [../]
[]

[Outputs]
  file_base = w_ro_0_out
  csv = true
  exodus = true
  [./console]
    type = Console
    perf_log = true
    output_linear = true
    max_rows = 25
  [../]

[]
