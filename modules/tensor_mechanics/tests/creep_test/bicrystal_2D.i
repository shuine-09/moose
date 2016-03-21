[Mesh]
  type = GeneratedMesh
  nx = 30
  ny = 15
  nz = 0
  xmin = 0
  xmax = 150
  ymin = 0
  ymax = 75
  zmin = 0
  zmax = 0
  dim = 2
  elem_type = QUAD4
  displacements = 'disp_x disp_y'
[]

[Variables]
  [./disp_x]
    block = 0
  [../]
  [./disp_y]
    block = 0
  [../]
[]

[GlobalParams]
  burgers_length  = 0.248e-3
  boltz_const = 1.38e-17
  temp = 800.0
[]

[AuxVariables]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./fp_yy]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./e_yy]
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
  [./glide_rate]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./athermal_resistance]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./euler1]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./euler2]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./euler3]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
[]

[Kernels]
  [./TensorMechanics]
    displacements = 'disp_x disp_y'
    use_displaced_mesh = true
  [../]
[]

[AuxKernels]
  [./stress_yy]
    type = RankTwoAux
    variable = stress_yy
    rank_two_tensor = stress
    index_j = 1
    index_i = 1
    execute_on = timestep_end
    block = 0
  [../]
  [./fp_yy]
    type = RankTwoAux
    variable = fp_yy
    rank_two_tensor = fp
    index_j = 1
    index_i = 1
    execute_on = timestep_end
    block = 0
  [../]
  [./e_yy]
    type = RankTwoAux
    variable = e_yy
    rank_two_tensor = lage
    index_j = 1
    index_i = 1
    execute_on = timestep_end
    block = 0
  [../]
  [./rho_m]
    type = MaterialStdVectorAux
    variable = rho_m
    property = mobile_dislocation_density
    index = 0
    execute_on = timestep_end
    block = 0
  [../]
  [./rho_i]
    type = MaterialStdVectorAux
    variable = rho_i
    property = immobile_dislocation_density
    index = 0
    execute_on = timestep_end
    block = 0
  [../]
  [./glide_rate]
    type = MaterialStdVectorAux
    variable = glide_rate
    property = glide_rate
    index = 0
    execute_on = timestep_end
    block = 0
  [../]
  [./athermal_slip_resistance]
    type = MaterialStdVectorAux
    variable = athermal_resistance
    property = athermal_slip_resistance
    index = 0
    execute_on = timestep_end
    block = 0
  [../]
  [./euler1]
    type = MaterialRealVectorValueAux
    variable = euler1
    property = Euler_angles
    component = 0
    execute_on = timestep_end
    block = 0
  [../]
  [./euler2]
    type = MaterialRealVectorValueAux
    variable = euler2
    property = Euler_angles
    component = 1
    execute_on = timestep_end
    block = 0
  [../]
  [./euler3]
    type = MaterialRealVectorValueAux
    variable = euler3
    property = Euler_angles
    component = 2
    execute_on = timestep_end
    block = 0
  [../]
[]

[BCs]
  [./symmy]
    type = PresetBC
    variable = disp_y
    boundary = bottom
    value = 0
  [../]
  [./symmx]
    type = PresetBC
    variable = disp_x
    boundary = left
    value = 0
  [../]
  [./dispy]
    type = FunctionPresetBC
    variable = disp_y
    boundary = top
    function = '75*0.0001*t'
  [../]
[]

[UserObjects]
  [./prop_read]
    type = ReadGrainCenterFromFile
    prop_file_name = 'euler_ang.txt'
    grain_center_file_name = 'grcenter_2D.txt'
    nprop = 3
    read_type = grain
    ngrain = 2
  [../]
  [./glide_rate]
    type = CPDislocationBasedGlideSlipRate
    variable_size = 12
    penalty_param = 1e-5
    disloc_line_dist = 1.5e-3
    jump_freq = 1e10
    active_enthal  = 6e-13
    exponentp = 0.28
    exponentq = 1.34
    prefactor = 1.0
    slip_sys_file_name = input_slip_sys.txt
    uo_mobile_dislocation_density_name = mobile_dislocation_density
    uo_thermal_slip_resistance_name = thermal_slip_resistance
    uo_athermal_slip_resistance_name = athermal_slip_resistance
  [../]
  [./thermal_slip_resistance]
    type = CPDislocationBasedThermalSlipResistance
    variable_size = 12
    thermal_resist = 500e-6
  [../]
  [./athermal_slip_resistance]
    type = CPDislocationBasedAthermalSlipResistance
    variable_size = 12
    rho_m_barrier_factor = 0.5
    shear_mod = 79731e-6
    rho_self_hard_factor = 1.0
    rho_latent_hard_factor = 0.2
    uo_mobile_dislocation_density_name = mobile_dislocation_density
    uo_immobile_dislocation_density_name = immobile_dislocation_density
  [../]
  [./mobile_dislocation_density]
    type = CrystalPlasticityStateVariable
    variable_size = 12
    groups = '0 12'  
    group_values = 20.0
    uo_state_var_evol_rate_comp_name = 'mobile_glide_rate_comp'
    scale_factor = '1.0'
  [../]
  [./immobile_dislocation_density]
    type = CrystalPlasticityStateVariable
    variable_size = 12
    groups = '0 12'  
    group_values = 30.0
    uo_state_var_evol_rate_comp_name = 'immobile_glide_rate_comp'
    scale_factor = '1.0'
  [../]
  [./mobile_glide_rate_comp]
    type = CPDislocationBasedMobileGlideRateComp
    variable_size = 12
    rho_mult_factor = 0.2
    rho_m_capture_radius = 1e-3
    rho_imm_factor = 0.05
    uo_mobile_dislocation_density_name = mobile_dislocation_density
    uo_immobile_dislocation_density_name = immobile_dislocation_density
    uo_glide_slip_rate_name = glide_rate
  [../]
  [./immobile_glide_rate_comp]
    type = CPDislocationBasedImmobileGlideRateComp
    variable_size = 12
    rho_imm_factor = 0.05
    uo_mobile_dislocation_density_name = mobile_dislocation_density
    uo_immobile_dislocation_density_name = immobile_dislocation_density
    uo_glide_slip_rate_name = glide_rate
  [../]
[]

[Materials]
  [./crysp]
    type = FiniteStrainUObasedCP
    block = 0
    stol = 1e-2
    rtol = 1e-6
    abs_tol = 1e-15
    tan_mod_type = none
    maximum_substep_iteration = 4
    uo_slip_rates = 'glide_rate'
    uo_slip_resistances = 'thermal_slip_resistance athermal_slip_resistance'
    uo_state_vars = 'mobile_dislocation_density immobile_dislocation_density'
    uo_state_var_evol_rate_comps = 'mobile_glide_rate_comp immobile_glide_rate_comp'
  [../]
  [./strain]
    type = ComputeFiniteStrain
    block = 0
    displacements = 'disp_x disp_y'
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensorCP
    block = 0
    C_ijkl = '0.279 0.12 0.12 0.279 0.12 0.279 0.08 0.08 0.08'
    fill_method = symmetric9
    read_prop_user_object = prop_read
  [../]
[]

[Postprocessors]
  [./stress_yy]
    type = ElementAverageValue
    variable = stress_yy
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./fp_yy]
    type = ElementAverageValue
    variable = fp_yy
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./e_yy]
    type = ElementAverageValue
    variable = e_yy
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
  [./glide_rate_aux]
    type = ElementAverageValue
    variable = glide_rate
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./athermal_resistance_aux]
    type = ElementAverageValue
    variable = athermal_resistance
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

  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 100'

#  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -ksp_gmres_restart -snes_ksp_ew_rtol0 -snes_ksp_ew_rtolmax -snes_ksp_ew_gamma -snes_ksp_ew_alpha -snes_ksp_ew_alpha2 -snes_ksp_ew_threshold'
#  petsc_options_value = 'lu superlu_dist 51 0.5 0.9 1 2 2 0.1'

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-10
  l_max_its = 60
  nl_max_its = 10

  dt = 1.0
#  dtmax = 10.0
#  dtmin = 1e-3

#  num_steps = 1
  end_time = 500.0
[]

[Outputs]
  csv = true
  exodus = true
[]

[Problem]
  use_legacy_uo_initialization = false
[]
