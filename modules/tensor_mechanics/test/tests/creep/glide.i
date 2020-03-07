[Mesh]
  type = GeneratedMesh
  nx = 30
  ny = 30
  nz = 30
  xmin = 0
  xmax = 60
  ymin = 0
  ymax = 60
  zmin = 0
  zmax = 60
  dim = 3
  elem_type = HEX8
  displacements = 'disp_x disp_y disp_z'
[]

[Variables]
  [./disp_x]
    block = 0
  [../]
  [./disp_y]
    block = 0
  [../]
  [./disp_z]
    block = 0
  [../]
[]

[GlobalParams]
  burgers_length  = 0.248e-3
  boltz_const = 1.38e-17
[]

[AuxVariables]
  [./stress_zz]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./strain_zz]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./fp_zz]
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
    displacements = 'disp_x disp_y disp_z'
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
  [./strain_zz]
    type = RankTwoAux
    variable = strain_zz
    rank_two_tensor = total_strain
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
  [./rho_m]
    type = MaterialStdVectorAux
    variable = rho_m
    property = mobile_dislocation_density
    index = 9
    execute_on = timestep_end
    block = 0
  [../]
  [./rho_i]
    type = MaterialStdVectorAux
    variable = rho_i
    property = immobile_dislocation_density
    index = 9
    execute_on = timestep_end
    block = 0
  [../]
  [./glide_rate]
    type = MaterialStdVectorAux
    variable = glide_rate
    property = glide_rate
    index = 9
    execute_on = timestep_end
    block = 0
  [../]
  [./athermal_slip_resistance]
    type = MaterialStdVectorAux
    variable = athermal_resistance
    property = athermal_slip_resistance
    index = 9
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
  [./symmz]
    type = PresetBC
    variable = disp_z
    boundary = back
    value = 0
  [../]
  #[./Pressure]
  #  [./load]
  #    #Applies the pressure
  #    boundary = front
  #    factor = 1.0
  #    function = -16e-6
  #    disp_x = disp_x
  #    disp_y = disp_y
  #    disp_z = disp_z
  #  [../]
  #[../]
  [./pull]
    type = FunctionPresetBC
    variable = disp_z
    boundary = front
    function = '0.2*t'
  [../]
[]

[UserObjects]
  [./prop_read]
    type = ElementPropertyReadFile
    prop_file_name = 'euler_angles.txt'
    nprop = 3
    read_type = voxel
  [../]
  [./glide_rate]
    type = CPDislocationBasedGlideSlipRate
    variable_size = 12
    penalty_param = 1e-4
    disloc_line_dist = 1.5e-3 #1.5e-3
    jump_freq = 1e10
    active_enthal  = 5.148e-13
    exponentp = 0.7
    exponentq = 1.4 #1.15
    prefactor = 1.0
    slip_incr_tol = 10.0
    slip_sys_file_name = input_slip_sys.txt
    uo_mobile_dislocation_density_name = mobile_dislocation_density
    uo_thermal_slip_resistance_name = thermal_slip_resistance
    uo_athermal_slip_resistance_name = athermal_slip_resistance
    #temp = 300.15 #
    temp = 573.15 #300C
    #temp = 773.15 #500C
    #temp = 973.15 #700C
    #temp = 1073.15 #800C
  [../]
  [./thermal_slip_resistance]
    type = CPDislocationBasedThermalSlipResistance
    variable_size = 12
    thermal_resist = 125e-6  # 514.8e-7
  [../]
  [./athermal_slip_resistance]
    type = CPDislocationBasedAthermalSlipResistance
    variable_size = 12
    rho_m_barrier_factor = 1.5
    shear_mod = 0.061 #0.0584
    rho_self_hard_factor = 1.2 #1.0
    rho_latent_hard_factor = 0.4 #0.2
    uo_mobile_dislocation_density_name = mobile_dislocation_density
    uo_immobile_dislocation_density_name = immobile_dislocation_density
  [../]
  [./mobile_dislocation_density]
    type = CrystalPlasticityStateVariable
    variable_size = 12
    groups = '0 12'
    group_values = 1
    uo_state_var_evol_rate_comp_name = 'mobile_glide_rate_comp'
    scale_factor = '1.0'
  [../]
  [./immobile_dislocation_density]
    type = CrystalPlasticityStateVariable
    variable_size = 12
    groups = '0 12'
    group_values = 1
    uo_state_var_evol_rate_comp_name = 'immobile_glide_rate_comp'
    scale_factor = '1.0'
  [../]
  [./mobile_glide_rate_comp]
    type = CPDislocationBasedMobileGlideRateComp
    variable_size = 12
    rho_mult_factor = 0.01 #0.143
    rho_m_capture_radius = 1e-3
    rho_imm_factor = 0.01 #0.025
    uo_mobile_dislocation_density_name = mobile_dislocation_density
    uo_immobile_dislocation_density_name = immobile_dislocation_density
    uo_glide_slip_rate_name = glide_rate
  [../]
  [./immobile_glide_rate_comp]
    type = CPDislocationBasedImmobileGlideRateComp
    variable_size = 12
    rho_imm_factor = 0.01
    uo_mobile_dislocation_density_name = mobile_dislocation_density
    uo_immobile_dislocation_density_name = immobile_dislocation_density
    uo_glide_slip_rate_name = glide_rate
  [../]
[]

[Materials]
  [./crysp]
    type = FiniteStrainUObasedCP
    block = 0
    stol = 1e-4
    rtol = 1e-6
    abs_tol = 1e-13
    tan_mod_type = exact
    zero_tol = 1.0e-13
    maximum_substep_iteration = 4
    maxiter = 100
    maxiter_state_variable = 100
    uo_slip_rates = 'glide_rate'
    uo_slip_resistances = 'thermal_slip_resistance athermal_slip_resistance'
    uo_state_vars = 'mobile_dislocation_density immobile_dislocation_density'
    uo_state_var_evol_rate_comps = 'mobile_glide_rate_comp immobile_glide_rate_comp'
  [../]
  [./strain]
    type = ComputeFiniteStrain
    block = 0
    displacements = 'disp_x disp_y disp_z'
    volumetric_locking_correction = false
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensorCP
    block = 0
    #C_ijkl = '0.1280 0.0812 0.0812 0.1280 0.0812 0.1280 0.0584 0.0584 0.0584'
    #C_ijkl = '0.284 0.121 0.121 0.284 0.121 0.284 0.081 0.081 0.081' #25C
    #C_ijkl = '0.261 0.111 0.111 0.261 0.111 0.111 0.075 0.075 0.075' #300C
    #C_ijkl = '0.243 0.104 0.104 0.243 0.104 0.104 0.07 0.07 0.07' #500C
    #C_ijkl = '0.223 0.095 0.095 0.223 0.095 0.095 0.064 0.064 0.064' #700C
    C_ijkl = '0.211 0.09 0.09 0.211 0.09 0.09 0.061 0.061 0.061' #800C
    fill_method = symmetric9
    read_prop_user_object = prop_read
  [../]
[]

[Postprocessors]
  [./stress_zz]
    type = ElementAverageValue
    variable = stress_zz
    block = 0
  [../]
  [./strain_zz]
    type = ElementAverageValue
    variable = strain_zz
    block = 0
  [../]
  [./fp_zz]
    type = ElementAverageValue
    variable = fp_zz
    block = 0
  [../]
  [./e_zz]
    type = ElementAverageValue
    variable = e_zz
    block = 0
  [../]
  [./rho_m]
    type = ElementAverageValue
    variable = rho_m
    block = 0
  [../]
  [./rho_i]
    type = ElementAverageValue
    variable = rho_i
    block = 0
  [../]
  [./glide_rate_aux]
    type = ElementAverageValue
    variable = glide_rate
    block = 0
  [../]
  [./athermal_resistance_aux]
    type = ElementAverageValue
    variable = athermal_resistance
    block = 0
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

  #Preconditioned JFNK (default#)
  solve_type = 'PJFNK'
#  solve_type = NEWTON

#  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
#  petsc_options_value = 'hypre boomeramg 100'

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu superlu_dist'

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6
  l_max_its = 20
  nl_max_its = 15

  #[./TimeStepper]
  #  type = FunctionDT
  #  time_dt = '0.1  0.1  1.0 1.0 10 10'
  #  time_t =  '0  1  11 100 101 1000'
  #[../]

  dt = 0.5
#  dtmax = 10.0
  dtmin = 0.01

#  num_steps = 1
  end_time = 50
[]

[Outputs]
  csv = true
  exodus = true
  #checkpoint = true
  file_base = T300_new_out2
[]

[Debug]
  show_var_residual_norms = true  ##Print residual at each iteration
[]
