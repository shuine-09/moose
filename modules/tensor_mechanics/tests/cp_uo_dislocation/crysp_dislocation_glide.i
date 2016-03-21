[Mesh]
  type = GeneratedMesh
  dim = 3
  elem_type = HEX8
  displacements = 'ux uy uz'
[]

[Variables]
  [./ux]
    block = 0
  [../]
  [./uy]
    block = 0
  [../]
  [./uz]
    block = 0
  [../]
[]

[AuxVariables]
  [./stress_zz]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./fp_zz]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./rotout]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./e_zz]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./rho_m]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./rho_i]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
[]

[Functions]
  [./tdisp]
    type = ParsedFunction
    value = 0.01*t
  [../]
[]

[Kernels]
  [./TensorMechanics]
    displacements = 'ux uy uz'
    use_displaced_mesh = true
  [../]
[]

[AuxKernels]
  [./stress_zz]
    type = RankTwoAux
    variable = stress_zz
    rank_two_tensor = stress
    index_j = 2
    index_i = 2
    execute_on = timestep_end
    block = 0
  [../]
  [./fp_zz]
    type = RankTwoAux
    variable = fp_zz
    rank_two_tensor = fp
    index_j = 2
    index_i = 2
    execute_on = timestep_end
    block = 0
  [../]
  [./e_zz]
    type = RankTwoAux
    variable = e_zz
    rank_two_tensor = lage
    index_j = 2
    index_i = 2
    execute_on = timestep_end
    block = 0
  [../]
  [./rotout]
    type = CrystalPlasticityRotationOutAux
    variable = rotout
    execute_on = timestep_end
    block = 0
  [../]
  [./rho_m]
    type = MaterialStdVectorAux
    variable = rho_m 
    property = mobile_dislocation_density
    index = 14
    execute_on = timestep_end
    block = 0
  [../]
  [./rho_i]
    type = MaterialStdVectorAux
    variable = rho_i 
    property = immobile_dislocation_density
    index = 14
    execute_on = timestep_end
    block = 0
  [../]
[]

[BCs]
  [./symmy]
    type = PresetBC
    variable = uy
    boundary = bottom
    value = 0
  [../]
  [./symmx]
    type = PresetBC
    variable = ux
    boundary = left
    value = 0
  [../]
  [./symmz]
    type = PresetBC
    variable = uz
    boundary = back
    value = 0
  [../]
  [./tdisp]
    type = FunctionPresetBC
    variable = uz
    boundary = front
    function = tdisp
  [../]
[]

[UserObjects]
  [./glide_slip_rate]
    type = CPDislocationBasedGlideSlipRate
    variable_size = 48
    penalty_param = 1e-3
    burgers_length  = 0.248e-6 #mm
    disloc_line_dist = 1.5e-6  #mm
    jump_freq = 1e6
    active_enthal  = 1.2161e-18
    temp = 323.0 
    slip_sys_file_name = input_slip_sys.txt
    uo_mobile_dislocation_density_name = mobile_dislocation_density
    uo_thermal_slip_resistance_name = thermal_slip_resistance
    uo_athermal_slip_resistance_name = athermal_slip_resistance
  [../]
  [./thermal_slip_resistance]
    type = CPDislocationBasedThermalSlipResistance
    variable_size = 48
    thermal_resist = 390.0
  [../]
  [./athermal_slip_resistance]
    type = CPDislocationBasedAthermalSlipResistance
    variable_size = 48
    burgers_length  = 0.248e-6 #mm
    uo_mobile_dislocation_density_name = mobile_dislocation_density
    uo_immobile_dislocation_density_name = immobile_dislocation_density
  [../]
  [./mobile_dislocation_density]
    type = CPDislocationBasedMobileDensity
    variable_size = 48
    initial_value = 2.0e7
    uo_state_var_evol_rate_comp_name = mobile_glide_rate_comp
    scale_factor = 1.0
  [../]
  [./immobile_dislocation_density]
    type = CPDislocationBasedImmobileDensity
    variable_size = 48
    initial_value = 3.0e7
    uo_state_var_evol_rate_comp_name = immobile_glide_rate_comp
    scale_factor = 1.0
  [../]
  [./mobile_glide_rate_comp]
    type = CPDislocationBasedMobileGlideRateComp
    variable_size = 48
    burgers_length  = 0.248e-6 #mm
    rho_imm_factor = 0.1
    uo_mobile_dislocation_density_name = mobile_dislocation_density
    uo_immobile_dislocation_density_name = immobile_dislocation_density
    uo_glide_slip_rate_name = glide_slip_rate
  [../]
  [./immobile_glide_rate_comp]
    type = CPDislocationBasedImmobileGlideRateComp
    variable_size = 48
    burgers_length  = 0.248e-6 #mm
    rho_imm_factor = 0.1
    uo_mobile_dislocation_density_name = mobile_dislocation_density
    uo_immobile_dislocation_density_name = immobile_dislocation_density
    uo_glide_slip_rate_name = glide_slip_rate
  [../]
[]

[Materials]
  [./crysp]
    type = FiniteStrainUObasedCP
    block = 0
    stol = 1.0e-5
    rtol = 1e-7
    abs_tol = 1e-15
    tan_mod_type = exact
    maximum_substep_iteration = 1
    uo_slip_rates = 'glide_slip_rate'
    uo_slip_resistances = 'thermal_slip_resistance athermal_slip_resistance'
    uo_state_vars = 'mobile_dislocation_density immobile_dislocation_density'
    uo_state_var_evol_rate_comps = 'mobile_glide_rate_comp immobile_glide_rate_comp'
  [../]
  [./strain]
    type = ComputeFiniteStrain
    block = 0
    displacements = 'ux uy uz'
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensorCP
    block = 0
    C_ijkl = '279.06e3 119.6e3 119.6e3 279.06e3 119.6e3 279.06e3 79.731e3 79.731e3 79.731e3'
    fill_method = symmetric9
  [../]
[]

[Postprocessors]
  [./stress_zz]
    type = ElementAverageValue
    variable = stress_zz
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./fp_zz]
    type = ElementAverageValue
    variable = fp_zz
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./e_zz]
    type = ElementAverageValue
    variable = e_zz
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./rho_m]
    type = ElementAverageValue
    variable = rho_m
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./rho_i]
    type = ElementAverageValue
    variable = rho_i
    block = 'ANY_BLOCK_ID 0'
    [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  petsc_options_iname = -pc_hypre_type
  petsc_options_value = boomerang
  #  nl_abs_tol = 1e-10
  nl_rel_step_tol = 1e-6
  dtmax = 10.0
  nl_rel_tol = 1e-6
  #  ss_check_tol = 1e-10
  end_time = 50.0
  dtmin = 1e-8
#  num_steps = 200
  dt = 0.1
  #  nl_abs_step_tol = 1e-10
  l_max_its = 60
  nl_max_its = 20
[]

[Outputs]
  file_base = out
  csv = true
  exodus = true
[]

[Problem]
  use_legacy_uo_initialization = false
[]
