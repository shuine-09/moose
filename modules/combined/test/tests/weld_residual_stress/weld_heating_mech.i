[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 120
  ny = 4
  xmax = .3
  ymax = .01
[]

[MeshModifiers]
  [./left_block]
    type = SubdomainBoundingBox
    block_id = 1
    bottom_left = '0 0 0'
    top_right = '.148 .01 0'
  [../]
  [./right_block]
    type = SubdomainBoundingBox
    block_id = 2
    bottom_left = '.152 0 0'
    top_right = '.3 .01 0'
  [../]
  [./bottom_weld]
    type = SubdomainBoundingBox
    block_id = 3
    bottom_left = '.148 0 0'
    top_right = '.152 .005 0'
  [../]
  [./top_weld]
    type = SubdomainBoundingBox
    block_id = 4
    bottom_left = '.148 .005 0'
    top_right = '.152 .01 0'
  [../]
  [./all_nodes]
    type = BoundingBoxNodeSet
    top_right = '.3 .01 0'
    bottom_left = '0 0 0'
    new_boundary = 'all_nodes'
  [../]
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
    use_displaced_mesh = true
    temperature = T
  [../]
  [./HeatSource1]
    type = WeldHeatSource
    variable = T
    block = 3
    value = 1e10
    weld_state = weld1
  [../]
  [./HeatSource2]
    type = WeldHeatSource
    variable = T
    block = 4
    value = 1e10
    weld_state = weld2
  [../]
[]

[BCs]
  [./left_fix_x]
    type = PresetBC
    variable = disp_x
    boundary = 'left right'
    value = 0.0
  [../]
  [./left_fix_y]
    type = PresetBC
    variable = disp_y
    boundary = 'left right'
    value = 0.0
  [../]
  [./left_temp]
    type = PresetBC
    variable = T
    boundary = 'left right'
    value = 20.0
  [../]
  # [./set_global_temp]
  #   type = DirichletBC
  #   variable = T
  #   boundary = 'all_nodes'
  #   value = '20.0'
  # [../]
  [./convective]    # Convective Start
    type = ConvectiveFluxFunction  # Convective flux, e.g. q'' = h*(Tw - Tf)
    boundary = 'left top bottom right'
    variable = T
    coefficient = 100
    T_infinity = 20
  [../]
[]

[Materials]
  [./non_weld_material]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 193e9
    poissons_ratio = 0.29
    block = '1 2'
  [../]
  [./non_weld_props]
    type = GenericConstantMaterial
    prop_names = 'thermal_conductivity specific_heat thermal_conductivity_dT'
    prop_values = '14 500 0'
    block = '1 2'
  [../]
  [./density]
    type = Density
    density = 799.0
  [../]
  [./thermal_strain]
    type = ComputeThermalExpansionEigenstrain
    stress_free_temperature = 20.0
    thermal_expansion_coeff = 1e-5
    eigenstrain_name = eigenstrain
    temperature = T
    block = '1 2'
  [../]

  [./weld_1]
    type = ComputeWeldingIsotropicElasticityTensor
    youngs_modulus = 193e9
    poissons_ratio = 0.29
    block = 3
    weld_state = weld1
  [../]
  [./hcm1]
    type = WeldingHeatConductionMaterial
    block = 3
    specific_heat = 500 #J/kg*C
    thermal_conductivity = 14 #W/m*C
    weld_state = weld1
  [../]
  [./thermal_strain1]
    type = ComputeWeldThermalExpansionEigenstrain
    stress_free_temperature = 1375.0
    thermal_expansion_coeff = 1e-5
    eigenstrain_name = eigenstrain
    temperature = T
    block = 3
    weld_state = weld1
  [../]

  [./weld_2]
    type = ComputeWeldingIsotropicElasticityTensor
    youngs_modulus = 193e9
    poissons_ratio = 0.29
    block = 4
    weld_state = weld2
  [../]
  [./hcm2]
    type = WeldingHeatConductionMaterial
    block = 4
    specific_heat = 500 #J/kg*C
    thermal_conductivity = 14 #W/m*C
    weld_state = weld2
  [../]
  [./thermal_strain2]
    type = ComputeWeldThermalExpansionEigenstrain
    stress_free_temperature = 1375.0
    thermal_expansion_coeff = 1e-5
    eigenstrain_name = eigenstrain
    temperature = T
    block = 4
    weld_state = weld2
  [../]

  [./strain]
    type = ComputeSmallStrain
    eigenstrain_names = eigenstrain
  [../]
  [./stress]
    type = ComputeLinearElasticStress
  [../]
[]

[UserObjects]
  [./weld1]
    type = WeldStateIndicator
    variable = T
    block = 3
    start_time = .1
    end_time = 1.0
    melting_temp = 1375 #C
    execute_on = 'timestep_end'
  [../]
  [./weld2]
    type = WeldStateIndicator
    variable = T
    block = 4
    start_time = 1.1
    end_time = 2.5
    melting_temp = 1375 #C
    execute_on = 'timestep_end'
  [../]
[]

# [Controls]
#   [./return_to_room_temp1]
#     type = TimePeriod
#     start_time = 1
#     end_time = 1.02
#     enable_objects = '*/set_global_temp'
#     execute_on = 'timestep_begin'
#     set_sync_times = true
#   [../]
# []

[Postprocessors]
  [./point_value1]
    type = PointValue
    point = '0.150 0.0 0.0'
    variable = T
  [../]
  [./point_value2]
    type = PointValue
    point = '0.150 0.01 0.0'
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
  scheme = bdf2
  solve_type = PJFNK

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu     superlu_dist'

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  l_tol = 1e-4
  dt = .01
  end_time = 4
[]

[Outputs]
  exodus = true
  csv = true
  print_perf_log = true
[]
