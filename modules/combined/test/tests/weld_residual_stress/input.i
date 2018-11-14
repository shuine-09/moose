[Mesh]
  type = FileMesh
  file = weld_mesh.e
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[AuxVariables]
  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./elastic_strain_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./elastic_strain_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./C1111]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./thermal_cond]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 0
  [../]
  [./stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
  [../]
 [./elastic_strain_xx]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = elastic_strain_xx
    index_i = 0
    index_j = 0
    execute_on = timestep_end
  [../]
  [./elastic_strain_yy]
     type = RankTwoAux
     rank_two_tensor = elastic_strain
     variable = elastic_strain_yy
     index_i = 1
     index_j = 1
     execute_on = timestep_end
   [../]
   [./C1111]
      type = RankFourAux
      rank_four_tensor = elasticity_tensor
      variable = C1111
      index_i = 0
      index_j = 0
      index_k = 0
      index_l = 0
      execute_on = timestep_end
   [../]
   [./thermal_cond]
     type = MaterialRealAux
     property = thermal_conductivity
     variable = thermal_cond
   [../]
[]

[Variables]
  [./T]
    initial_condition = 20.0
  [../]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
[]

[Kernels]
  [./HeatDiff]
    type = HeatConduction
    variable = T
  [../]
  [./HeatTdot]
    type = HeatConductionTimeDerivative
    variable = T
  [../]
  [./TensorMechanics]
    use_displaced_mesh = false
    temperature = T
  [../]
  [./HeatSource1]
    type = WeldHeatSource
    variable = T
    block = 1
    value = 1e11
    weld_state = weld1
  [../]
  [./HeatSource2]
    type = WeldHeatSource
    variable = T
    block = 2
    value = 1e11
    weld_state = weld2
  [../]
  [./HeatSource3]
    type = WeldHeatSource
    variable = T
    block = 3
    value = 1e11
    weld_state = weld3
  [../]
  [./HeatSource4]
    type = WeldHeatSource
    variable = T
    block = 4
    value = 1e11
    weld_state = weld4
  [../]
  [./HeatSource5]
    type = WeldHeatSource
    variable = T
    block = 5
    value = 1e11
    weld_state = weld5
  [../]
  [./HeatSource6]
    type = WeldHeatSource
    variable = T
    block = 6
    value = 1e11
    weld_state = weld6
  [../]
  [./HeatSource7]
    type = WeldHeatSource
    variable = T
    block = 7
    value = 1e11
    weld_state = weld7
  [../]
  [./HeatSource8]
    type = WeldHeatSource
    variable = T
    block = 8
    value = 1e11
    weld_state = weld8
  [../]
  [./HeatSource9]
    type = WeldHeatSource
    variable = T
    block = 9
    value = 1e11
    weld_state = weld9
  [../]
  [./HeatSource10]
    type = WeldHeatSource
    variable = T
    block = 10
    value = 1e11
    weld_state = weld10
  [../]
  [./HeatSource11]
    type = WeldHeatSource
    variable = T
    block = 11
    value = 1e11
    weld_state = weld11
  [../]
  [./HeatSource12]
    type = WeldHeatSource
    variable = T
    block = 12
    value = 1e11
    weld_state = weld12
  [../]
[]

[BCs]
  [./left_fix_x]
    type = PresetBC
    variable = disp_x
    boundary = '1 3 4'
    value = 0.0
  [../]
  [./left_fix_y]
    type = PresetBC
    variable = disp_y
    boundary = '1'
    value = 0.0
  [../]
  # [./set_global_temp]
  #   type = DirichletBC
  #   variable = T
  #   boundary = '1000'
  #   value = '20.0'
  # [../]
[]

[Materials]
  ### Non-Weld ###
  [./non_weld_material]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 200e9
    poissons_ratio = 0.29
    block = '100'
  [../]
  [./non_weld_props]
    type = GenericConstantMaterial
    prop_names = 'thermal_conductivity specific_heat thermal_conductivity_dT'
    prop_values = '14.6 450 0'
    block = '100'
  [../]

  ### Welds ###
  [./weld_1]
    type = ComputeWeldingIsotropicElasticityTensor
    youngs_modulus = 200e9
    poissons_ratio = 0.29
    block = '1'
    weld_state = weld1
  [../]
  [./weld_2]
    type = ComputeWeldingIsotropicElasticityTensor
    youngs_modulus = 200e9
    poissons_ratio = 0.29
    block = '2'
    weld_state = weld2
  [../]
  [./weld_3]
    type = ComputeWeldingIsotropicElasticityTensor
    youngs_modulus = 200e9
    poissons_ratio = 0.29
    block = '3'
    weld_state = weld3
  [../]
  [./weld_4]
    type = ComputeWeldingIsotropicElasticityTensor
    youngs_modulus = 200e9
    poissons_ratio = 0.29
    block = '4'
    weld_state = weld4
  [../]
  [./weld_5]
    type = ComputeWeldingIsotropicElasticityTensor
    youngs_modulus = 200e9
    poissons_ratio = 0.29
    block = '5'
    weld_state = weld5
  [../]
  [./weld_6]
    type = ComputeWeldingIsotropicElasticityTensor
    youngs_modulus = 200e9
    poissons_ratio = 0.29
    block = '6'
    weld_state = weld6
  [../]
  [./weld_7]
    type = ComputeWeldingIsotropicElasticityTensor
    youngs_modulus = 200e9
    poissons_ratio = 0.29
    block = '7'
    weld_state = weld7
  [../]
  [./weld_8]
    type = ComputeWeldingIsotropicElasticityTensor
    youngs_modulus = 200e9
    poissons_ratio = 0.29
    block = '8'
    weld_state = weld8
  [../]
  [./weld_9]
    type = ComputeWeldingIsotropicElasticityTensor
    youngs_modulus = 200e9
    poissons_ratio = 0.29
    block = '9'
    weld_state = weld9
  [../]
  [./weld_10]
    type = ComputeWeldingIsotropicElasticityTensor
    youngs_modulus = 200e9
    poissons_ratio = 0.29
    block = '10'
    weld_state = weld10
  [../]
  [./weld_11]
    type = ComputeWeldingIsotropicElasticityTensor
    youngs_modulus = 200e9
    poissons_ratio = 0.29
    block = '11'
    weld_state = weld11
  [../]
  [./weld_12]
    type = ComputeWeldingIsotropicElasticityTensor
    youngs_modulus = 200e9
    poissons_ratio = 0.29
    block = '12'
    weld_state = weld12
  [../]

  [./hcm_1]
    type = WeldingHeatConductionMaterial
    block = '1'
    specific_heat = 450 #J/kg*K
    thermal_conductivity = 14.6 #W/m*K
    weld_state = weld1
  [../]
  [./hcm_2]
    type = WeldingHeatConductionMaterial
    block = '2'
    specific_heat = 450 #J/kg*K
    thermal_conductivity = 14.6 #W/m*K
    weld_state = weld2
  [../]
  [./hcm_3]
    type = WeldingHeatConductionMaterial
    block = '3'
    specific_heat = 450 #J/kg*K
    thermal_conductivity = 14.6 #W/m*K
    weld_state = weld3
  [../]
  [./hcm_4]
    type = WeldingHeatConductionMaterial
    block = '4'
    specific_heat = 450 #J/kg*K
    thermal_conductivity = 14.6 #W/m*K
    weld_state = weld4
  [../]
  [./hcm_5]
    type = WeldingHeatConductionMaterial
    block = '5'
    specific_heat = 450 #J/kg*K
    thermal_conductivity = 14.6 #W/m*K
    weld_state = weld5
  [../]
  [./hcm_6]
    type = WeldingHeatConductionMaterial
    block = '6'
    specific_heat = 450 #J/kg*K
    thermal_conductivity = 14.6 #W/m*K
    weld_state = weld6
  [../]
  [./hcm_7]
    type = WeldingHeatConductionMaterial
    block = '7'
    specific_heat = 450 #J/kg*K
    thermal_conductivity = 14.6 #W/m*K
    weld_state = weld7
  [../]
  [./hcm_8]
    type = WeldingHeatConductionMaterial
    block = '8'
    specific_heat = 450 #J/kg*K
    thermal_conductivity = 14.6 #W/m*K
    weld_state = weld8
  [../]
  [./hcm_9]
    type = WeldingHeatConductionMaterial
    block = '9'
    specific_heat = 450 #J/kg*K
    thermal_conductivity = 14.6 #W/m*K
    weld_state = weld9
  [../]
  [./hcm_10]
    type = WeldingHeatConductionMaterial
    block = '10'
    specific_heat = 450 #J/kg*K
    thermal_conductivity = 14.6 #W/m*K
    weld_state = weld10
  [../]
  [./hcm_11]
    type = WeldingHeatConductionMaterial
    block = '11'
    specific_heat = 450 #J/kg*K
    thermal_conductivity = 14.6 #W/m*K
    weld_state = weld11
  [../]
  [./hcm_12]
    type = WeldingHeatConductionMaterial
    block = '12'
    specific_heat = 450 #J/kg*K
    thermal_conductivity = 14.6 #W/m*K
    weld_state = weld12
  [../]

  [./thermal_strain1]
    type = ComputeWeldThermalExpansionEigenstrain
    stress_free_temperature = 1673
    thermal_expansion_coeff = 19.5e-6
    eigenstrain_name = eigenstrain
    temperature = T
    block = 1
    weld_state = weld1
  [../]
  [./thermal_strain2]
    type = ComputeWeldThermalExpansionEigenstrain
    stress_free_temperature = 1673
    thermal_expansion_coeff = 19.5e-6
    eigenstrain_name = eigenstrain
    temperature = T
    block = 2
    weld_state = weld2
  [../]
  [./thermal_strain3]
    type = ComputeWeldThermalExpansionEigenstrain
    stress_free_temperature = 1673
    thermal_expansion_coeff = 19.5e-6
    eigenstrain_name = eigenstrain
    temperature = T
    block = 3
    weld_state = weld3
  [../]
  [./thermal_strain4]
    type = ComputeWeldThermalExpansionEigenstrain
    stress_free_temperature = 1673
    thermal_expansion_coeff = 19.5e-6
    eigenstrain_name = eigenstrain
    temperature = T
    block = 4
    weld_state = weld4
  [../]
  [./thermal_strain5]
    type = ComputeWeldThermalExpansionEigenstrain
    stress_free_temperature = 1673
    thermal_expansion_coeff = 19.5e-6
    eigenstrain_name = eigenstrain
    temperature = T
    block = 5
    weld_state = weld5
  [../]
  [./thermal_strain6]
    type = ComputeWeldThermalExpansionEigenstrain
    stress_free_temperature = 1673
    thermal_expansion_coeff = 19.5e-6
    eigenstrain_name = eigenstrain
    temperature = T
    block = 6
    weld_state = weld6
  [../]
  [./thermal_strain7]
    type = ComputeWeldThermalExpansionEigenstrain
    stress_free_temperature = 1673
    thermal_expansion_coeff = 19.5e-6
    eigenstrain_name = eigenstrain
    temperature = T
    block = 7
    weld_state = weld7
  [../]
  [./thermal_strain8]
    type = ComputeWeldThermalExpansionEigenstrain
    stress_free_temperature = 1673
    thermal_expansion_coeff = 19.5e-6
    eigenstrain_name = eigenstrain
    temperature = T
    block = 8
    weld_state = weld8
  [../]
  [./thermal_strain9]
    type = ComputeWeldThermalExpansionEigenstrain
    stress_free_temperature = 1673
    thermal_expansion_coeff = 19.5e-6
    eigenstrain_name = eigenstrain
    temperature = T
    block = 9
    weld_state = weld9
  [../]
  [./thermal_strain10]
    type = ComputeWeldThermalExpansionEigenstrain
    stress_free_temperature = 1673
    thermal_expansion_coeff = 19.5e-6
    eigenstrain_name = eigenstrain
    temperature = T
    block = 10
    weld_state = weld10
  [../]
  [./thermal_strain11]
    type = ComputeWeldThermalExpansionEigenstrain
    stress_free_temperature = 1673
    thermal_expansion_coeff = 19.5e-6
    eigenstrain_name = eigenstrain
    temperature = T
    block = 11
    weld_state = weld11
  [../]
  [./thermal_strain12]
    type = ComputeWeldThermalExpansionEigenstrain
    stress_free_temperature = 1673
    thermal_expansion_coeff = 19.5e-6
    eigenstrain_name = eigenstrain
    temperature = T
    block = 12
    weld_state = weld12
  [../]

  ### Everything ####
  [./strain]
    type = ComputeIncrementalSmallStrain
    eigenstrain_names = eigenstrain
  [../]
  [./thermal_strain]
    type = ComputeThermalExpansionEigenstrain
    stress_free_temperature = 293 #1673C
    thermal_expansion_coeff = 19.5e-6
    # thermal_expansion_coeff = 1.0e-6
    eigenstrain_name = eigenstrain
    temperature = T
    block = 100
  [../]
  [./stress]
    type = ComputeFiniteStrainElasticStress
  [../]
  [./density]
    type = Density
    density = 8029.0
  [../]
[]

[UserObjects]
  [./weld1]
    type = WeldStateIndicator
    variable = T
    block = 1
    start_time = .01
    end_time = 2.0
    melting_temp = 1673 #C
    execute_on = 'timestep_begin'
  [../]
  [./weld2]
    type = WeldStateIndicator
    variable = T
    block = 2
    start_time = 2.01
    end_time = 4
    melting_temp = 1673 #C
    execute_on = 'timestep_begin'
  [../]
  [./weld3]
    type = WeldStateIndicator
    variable = T
    block = 3
    start_time = 4.01
    end_time = 6.0
    melting_temp = 1673 #C
    execute_on = 'timestep_begin'
  [../]
  [./weld4]
    type = WeldStateIndicator
    variable = T
    block = 4
    start_time = 6.01
    end_time = 8
    melting_temp = 1673 #C
    execute_on = 'timestep_begin'
  [../]
  [./weld5]
    type = WeldStateIndicator
    variable = T
    block = 5
    start_time = 10.01
    end_time = 12.0
    melting_temp = 1673 #C
    execute_on = 'timestep_begin'
  [../]
  [./weld6]
    type = WeldStateIndicator
    variable = T
    block = 6
    start_time = 12.01
    end_time = 14.0
    melting_temp = 1673 #C
    execute_on = 'timestep_begin'
  [../]
  [./weld7]
    type = WeldStateIndicator
    variable = T
    block = 7
    start_time = 14.01
    end_time = 16.0
    melting_temp = 1673 #C
    execute_on = 'timestep_begin'
  [../]
  [./weld8]
    type = WeldStateIndicator
    variable = T
    block = 8
    start_time = 16.01
    end_time = 18.0
    melting_temp = 1673 #C
    execute_on = 'timestep_begin'
  [../]
  [./weld9]
    type = WeldStateIndicator
    variable = T
    block = 9
    start_time = 18.01
    end_time = 20.0
    melting_temp = 1673 #C
    execute_on = 'timestep_begin'
  [../]
  [./weld10]
    type = WeldStateIndicator
    variable = T
    block = 10
    start_time =20.01
    end_time = 22.0
    melting_temp = 1673 #C
    execute_on = 'timestep_begin'
  [../]
  [./weld11]
    type = WeldStateIndicator
    variable = T
    block = 11
    start_time = 22.01
    end_time = 24.0
    melting_temp = 1673 #C
    execute_on = 'timestep_begin'
  [../]
  [./weld12]
    type = WeldStateIndicator
    variable = T
    block = 12
    start_time = 24.01
    end_time = 26.0
    melting_temp = 1673 #C
    execute_on = 'timestep_begin'
  [../]
[]

# [Controls]
#   [./return_to_room_temp1]
#     type = TimePeriod
#     start_time = 1.99
#     end_time = 2.0
#     enable_objects = '*/set_global_temp'
#     execute_on = 'timestep_begin'
#     set_sync_times = true
#   [../]
#   [./return_to_room_temp2]
#     type = TimePeriod
#     start_time = 3.99
#     end_time = 4.0
#     enable_objects = '*/set_global_temp'
#     execute_on = 'timestep_begin'
#     set_sync_times = true
#   [../]
#   [./return_to_room_temp3]
#     type = TimePeriod
#     start_time = 5.99
#     end_time = 6.0
#     enable_objects = '*/set_global_temp'
#     execute_on = 'timestep_begin'
#     set_sync_times = true
#   [../]
#   [./return_to_room_temp4]
#     type = TimePeriod
#     start_time = 7.99
#     end_time = 8.0
#     enable_objects = '*/set_global_temp'
#     execute_on = 'timestep_begin'
#     set_sync_times = true
#   [../]
#   [./return_to_room_temp5]
#     type = TimePeriod
#     start_time = 9.99
#     end_time = 10.0
#     enable_objects = '*/set_global_temp'
#     execute_on = 'timestep_begin'
#     set_sync_times = true
#   [../]
#   [./return_to_room_temp6]
#     type = TimePeriod
#     start_time = 11.99
#     end_time = 12.0
#     enable_objects = '*/set_global_temp'
#     execute_on = 'timestep_begin'
#     set_sync_times = true
#   [../]
#   [./return_to_room_temp7]
#     type = TimePeriod
#     start_time = 13.99
#     end_time = 14.0
#     enable_objects = '*/set_global_temp'
#     execute_on = 'timestep_begin'
#     set_sync_times = true
#   [../]
#   [./return_to_room_temp8]
#     type = TimePeriod
#     start_time = 15.99
#     end_time = 16.0
#     enable_objects = '*/set_global_temp'
#     execute_on = 'timestep_begin'
#     set_sync_times = true
#   [../]
#   [./return_to_room_temp9]
#     type = TimePeriod
#     start_time = 17.99
#     end_time = 18.0
#     enable_objects = '*/set_global_temp'
#     execute_on = 'timestep_begin'
#     set_sync_times = true
#   [../]
#   [./return_to_room_temp10]
#     type = TimePeriod
#     start_time = 19.99
#     end_time = 20.0
#     enable_objects = '*/set_global_temp'
#     execute_on = 'timestep_begin'
#     set_sync_times = true
#   [../]
#   [./return_to_room_temp11]
#     type = TimePeriod
#     start_time = 21.99
#     end_time = 22.0
#     enable_objects = '*/set_global_temp'
#     execute_on = 'timestep_begin'
#     set_sync_times = true
#   [../]
#   [./return_to_room_temp12]
#     type = TimePeriod
#     start_time = 23.99
#     end_time = 24.0
#     enable_objects = '*/set_global_temp'
#     execute_on = 'timestep_begin'
#     set_sync_times = true
#   [../]
#   [./return_to_room_temp_end]
#     type = TimePeriod
#     start_time = 25.99
#     end_time = 26.0
#     enable_objects = '*/set_global_temp'
#     execute_on = 'timestep_end'
#     set_sync_times = true
#   [../]
# []

[Postprocessors]
  [./point_1]
    type = PointValue
    point = '2.358405 0.830387 0.0'
    variable = T
  [../]
  [./point_2]
    type = PointValue
    point = '2.35143 0.833154 0.0'
    variable = T
  [../]
  [./point_3]
    type = PointValue
    point = '2.343521 0.834805 0.0'
    variable = T
  [../]
  [./point_4]
    type = PointValue
    point = '2.363986 0.82951 0.0'
    variable = T
  [../]
  [./point_5]
    type = PointValue
    point = '2.368988 0.831126 0.0'
    variable = T
  [../]
  [./point_6]
    type = PointValue
    point = '2.37399 0.832969 0.0'
    variable = T
  [../]
  [./point_7]
    type = PointValue
    point = '2.376928 0.831615 0.0'
    variable = T
  [../]
  [./point_8]
    type = PointValue
    point = '2.376928 0.835385 0.0'
    variable = T
  [../]
  [./point_9]
    type = PointValue
    point = '2.330617 0.830866 0.0'
    variable = T
  [../]
  [./point_10]
    type = PointValue
    point = '2.340599 0.842133 0.0'
    variable = T
  [../]
  [./point_11]
    type = PointValue
    point = '2.333676 0.844261 0.0'
    variable = T
  [../]
  [./point_12]
    type = PointValue
    point = '2.322176 0.832301 0.0'
    variable = T
  [../]
[]


[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  # scheme = bdf2
  solve_type = PJFNK

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu     superlu_dist'
  l_max_its = 100
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  l_tol = 1e-4
  dt = .01
  end_time = 26.0
[]

[Outputs]
  exodus = true
  csv = true
  print_perf_log = true
[]
