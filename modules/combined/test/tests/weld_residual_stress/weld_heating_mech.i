[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 20
  ny = 10
  xmax = 10
  ymax = 5
[]

[MeshModifiers]
  [./add_block1]
    type = SubdomainBoundingBox
     block_id = 1
     bottom_left = '0 0 0'
     top_right = '5 5 0'
  [../]
  [./add_block2]
    type = SubdomainBoundingBox
     block_id = 2
     bottom_left = '5 0 0'
     top_right = '10 5 0'
  [../]
[]



[GlobalParams]
  block = '1 2'
  displacements = 'disp_x disp_y'
[]

[AuxVariables]
  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./elastic_strain_xx]
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

 [./elastic_strain_xx]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = elastic_strain_xx
    index_i = 0
    index_j = 0
    execute_on = timestep_end
  [../]
[]

[Variables]
  [./T]
    initial_condition = 300.0
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
  [./HeatSource1]
    type = WeldHeatSource
    variable = T
    block = 1
    value = 100
    weld_state = weld1
  [../]
  [./HeatSource2]
    type = WeldHeatSource
    variable = T
    block = 2
    value = 100
    weld_state = weld2
  [../]
  [./TensorMechanics]
    use_displaced_mesh = true
    temperature = T
  [../]
[]

[BCs]
  [./lefttemp]
    type = DirichletBC
    boundary = left
    variable = T
    value = 400
  [../]
  [./righttemp]
    type = DirichletBC
    boundary = right
    variable = T
    value = 300
  [../]
  [./left_fix_x]
    type = PresetBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]
  [./left_fix_y]
    type = PresetBC
    variable = disp_y
    boundary = left
    value = 0.0
  [../]
  [./right_x]
    type = FunctionPresetBC
    variable = disp_x
    boundary = right
    function = '0.01 * t'
  [../]
[]

[Materials]
  [./k]
    type = GenericConstantMaterial
    prop_names = 'thermal_conductivity'
    prop_values = '0.95' #copper in cal/(cm sec C)
  [../]
  [./cp]
    type = GenericConstantMaterial
    prop_names = 'specific_heat'
    prop_values = '0.092' #copper in cal/(g C)
  [../]
  [./rho]
    type = GenericConstantMaterial
    prop_names = 'density'
    prop_values = '8.92' #copper in g/(cm^3)
  [../]

  [./elasticity_tensor1]
    type = ComputeWeldingIsotropicElasticityTensor
    youngs_modulus = 1.0e6
    poissons_ratio = 0.3
    block = 1
    weld_state = weld1
  [../]

  [./elasticity_tensor2]
    type = ComputeWeldingIsotropicElasticityTensor
    youngs_modulus = 1.0e6
    poissons_ratio = 0.3
    block = 2
    weld_state = weld2
  [../]

  [./strain]
    type = ComputeSmallStrain
    eigenstrain_names = eigenstrain
  [../]
  [./thermal_strain]
    type = ComputeThermalExpansionEigenstrain
    stress_free_temperature = 0.0
    thermal_expansion_coeff = 1e-5
    eigenstrain_name = eigenstrain
    temperature = T
  [../]
  [./stress]
    type = ComputeLinearElasticStress
  [../]
[]

[UserObjects]
  [./weld1]
    type = WeldStateIndicator
    variable = T
    block = 1
    start_time = 1
    end_time = 100
    melting_temp = 400
    execute_on = 'timestep_begin'
  [../]
  [./weld2]
    type = WeldStateIndicator
    variable = T
    block = 2
    start_time = 3
    end_time = 100
    melting_temp = 400
    execute_on = 'timestep_begin'
  [../]
[]

[Postprocessors]
  [./point_value1]
    type = PointValue
    point = '2.5 2.5 0'
    variable = T
  [../]
  [./point_value2]
    type = PointValue
    point = '7.5 2.5 0'
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

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  l_tol = 1e-4
  dt = 1
  end_time = 5
[]

[Outputs]
  exodus = true
  print_perf_log = true
[]
