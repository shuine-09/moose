[Mesh]
  type = MeshGeneratorMesh
  uniform_refine = 0
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[MeshGenerators]
  [./mesh]
    type = FileMeshGenerator
    file = 2d.e
  [../]
  [./add_set1]
    type = SubdomainBoundingBoxGenerator
    input = mesh
    block_id = 2
    bottom_left = '0 0 0'
    top_right = '100 10 0'
  [../]
  [./add_set2]
    type = SubdomainBoundingBoxGenerator
    input = add_set1
    block_id = 1
    bottom_left = '0 10 0'
    top_right = '100 40 0'
  [../]
  [./add_bnd]
    type = SideSetsFromAllElementFaces
    input = add_set2
    block = 1
    new_boundary = 'moving_boundary'
  [../]
[]

[Variables]
  # [./disp_x]
  # [../]
  # [./disp_y]
  # [../]
  [./temp]
    initial_condition = 300
  [../]
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
  [./temp_aux]
  [../]
  [./activated_elem]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  [./time]
    type = ADHeatConductionTimeDerivative
    variable = temp
  [../]
  [./heat_conduct]
    type = ADHeatConduction
    variable = temp
    use_displaced_mesh = true
  [../]
  [./heat_source]
    type = ADMatHeatSource
    material_property = volumetric_heat
    variable = temp
    scalar = 1
    use_displaced_mesh = true
  [../]
  [./disp_x]
    type = ADStressDivergenceTensors
    component = 0
    variable = disp_x
    use_displaced_mesh = true
  [../]
  [./disp_y]
    type = ADStressDivergenceTensors
    component = 1
    variable = disp_y
    use_displaced_mesh = true
  [../]
[]


[AuxKernels]
  [./activated_elem]
    type = ActivatedElementsMarker
    melt_temperature = 700
    temp_aux = temp_aux
    variable = activated_elem
    execute_on = timestep_begin
    #marker_uo = activated_elem_uo
  [../]
  [./stress_xx]
    type = ADRankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 0
  [../]

 [./elastic_strain_xx]
    type = ADRankTwoAux
    rank_two_tensor = elastic_strain
    variable = elastic_strain_xx
    index_i = 0
    index_j = 0
    execute_on = timestep_end
  [../]
[]

[UserObjects]
  [./activated_elem_uo]
    type = ActivatedElementsMarkerUO
    melt_temperature = 700
    temp_aux = temp_aux
    execute_on = timestep_begin
  [../]
[]

[BCs]
  [./bottom_fix_x]
    type = PresetBC
    variable = disp_x
    boundary = 3
    value = 0.0
  [../]
  [./bottom_fix_y]
    type = PresetBC
    variable = disp_y
    boundary = 3
    value = 0.0
  [../]
  # [./fix_y]
  #   type = PresetBC
  #   variable = disp_y
  #   boundary = '101 102'
  #   value = 0.0
  # [../]
  # [./fix_x]
  #   type = PresetBC
  #   variable = disp_y
  #   boundary = 101
  #   value = 0.0
  # [../]

  [./convective]
    type = ConvectiveFluxFunction # Convective flux, e.g. q'' = h*(Tw - Tf)
    variable = temp
    boundary = 'moving_boundary'
    coefficient = 1e-2
    T_infinity = 300
    marker_uo = activated_elem_uo
  [../]
  [./convective_substrate]
    type = ConvectiveFluxFunction # Convective flux, e.g. q'' = h*(Tw - Tf)
    variable = temp
    boundary = '1 2 3 4 5'
    coefficient = 1e-2
    T_infinity = 300
  [../]
[]

[Functions]
  [./yield2]
    type = PiecewiseLinear
    x = '400 500 600'
    y = '6e8 5e8 4e8'
  [../]
  [./yield1]
    type = PiecewiseLinear
    x = '400 500 600'
    y = '6e4 6e4 6e4'
  [../]
[]

[Modules/TensorMechanics/Master]
  [./all]
    strain = FINITE
    incremental = true
    add_variables = true
    generate_output = 'stress_yy plastic_strain_xx plastic_strain_yy'
    use_automatic_differentiation = true
    eigenstrain_names = eigenstrain
  [../]
[]

[Materials]
  [./youngs_modulus]
    type = ADPiecewiseLinearInterpolationMaterial
    xy_data = '0          10e+8
               999.9999   10e+8
               1000       9.94e+6
               99900      10e5'
    property = youngs_modulus
    variable = temp
  [../]

  [./elasticity_tensor]
    type = ADComputeVariableIsotropicElasticityTensor
    youngs_modulus = youngs_modulus
    #youngs_modulus = 2e8
    poissons_ratio = 0.3
  [../]
  # [./strain]
  #   type = ADComputeFiniteStrain
  #   eigenstrain_names = eigenstrain
  # [../]
  # [./stress]
  #   type = ADComputeFiniteStrainElasticStress
  # [../]
  [./plas1]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plasticity1'
    max_iterations = 100
    absolute_tolerance = 1e-8
    relative_tolerance = 1e-8
    block = 1
  [../]
  [./plasticity1]
    type = ADIsotropicPlasticityStressUpdate
    hardening_constant = 0
    yield_stress_function = yield1
    temperature = temp
    block = 1
  [../]
  [./plas2]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'plasticity2'
    max_iterations = 100
    absolute_tolerance = 1e-8
    relative_tolerance = 1e-8
    block = 2
  [../]
  [./plasticity2]
    type = ADIsotropicPlasticityStressUpdate
    hardening_constant = 0
    yield_stress_function = yield2
    temperature = temp
    block = 2
  [../]
  [./thermal]
    type = ADComputeThermalExpansionEigenstrain
    stress_free_temperature = 300
    thermal_expansion_coeff = 8e-8
    temperature = temp
    eigenstrain_name = eigenstrain
    activated_elem_aux = activated_elem
  [../]
  [./heat]
    type = ADHeatConductionMaterial
    specific_heat = 603
    thermal_conductivity = 0.5
  [../]
  [./volumetric_heat]
    type = DoubleEllipsoidHeatSource
    a = 6
    b = 3
    c = 1.5
    power = 600
    efficienty = 0.45
    factor = 1
    velocity = 10
  [../]
  [./density]
    type = ADDensity
    density = 4.43e-6
  [../]
[]

[Preconditioning]
  [./full]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient

  solve_type = NEWTON

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6
  l_tol = 1e-3
  line_search = 'none'
  nl_max_its = 15

  petsc_options_iname = '-ksp_type -pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'preonly lu       superlu_dist'

  l_max_its = 100
[]

[Postprocessors]
  [./elastic_strain_xx]
    type = ElementAverageValue
    variable = elastic_strain_xx
  [../]
  [./elastic_stress_xx]
    type = ElementAverageValue
    variable = stress_xx
  [../]
  [./temp]
    type = AverageNodalVariableValue
    variable = temp
  [../]
[]

[Outputs]
  exodus = true
[]
