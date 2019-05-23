[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  type = MeshGeneratorMesh
[]

[MeshGenerators]
  [./mesh]
    type = FileMeshGenerator
    file = hollow.e
  [../]
  [./add_set1]
    type = SubdomainBoundingBoxGenerator
    input = mesh
    block_id = 2
    bottom_left = '-50 -50 0'
    top_right = '50 50 5'
  [../]
  [./add_set2]
    type = SubdomainBoundingBoxGenerator
    input = add_set1
    block_id = 3
    bottom_left = '50 50 5'
    top_right = '50 50 20'
  [../]
  [./add_bnd]
    type = SideSetsFromAllElementFaces
    input = add_set2
    block = 1
    new_boundary = 'moving_boundary'
  [../]
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./disp_z]
  [../]
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
  [stress_radial]
    family = MONOMIAL
    order = CONSTANT
  []
  [stress_hoop]
    family = MONOMIAL
    order = CONSTANT
  []
  [stress_axial]
    family = MONOMIAL
    order = CONSTANT
  []
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
  [./disp_z]
    type = ADStressDivergenceTensors
    component = 2
    variable = disp_z
    use_displaced_mesh = true
  [../]
[]


[AuxKernels]
  [./activated_elem]
    type = ActivatedElementsMarker
    melt_temperature = 480
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
  [./radial]
    type = ADRankTwoScalarAux
    variable = stress_radial
    rank_two_tensor = stress
    scalar_type = RadialStress
    point1 = '0 0 5'
    point2 = '0 1 5'
    execute_on = timestep_end
  [../]
  [./hoop]
    type = ADRankTwoScalarAux
    variable = stress_hoop
    rank_two_tensor = stress
    scalar_type = HoopStress
    point1 = '0 0 5'
    point2 = '0 1 5'
    execute_on = timestep_end
  [../]
  [./axial]
    type = ADRankTwoScalarAux
    variable = stress_axial
    rank_two_tensor = stress
    scalar_type = AxialStress
    point1 = '0 0 5'
    point2 = '0 1 5'
    execute_on = timestep_end
  [../]
[]

[UserObjects]
  [./activated_elem_uo]
    type = ActivatedElementsMarkerUO
    melt_temperature = 480
    temp_aux = temp_aux
    execute_on = timestep_begin
  [../]
[]

[BCs]
  [./bottom_fix_x]
    type = PresetBC
    variable = disp_x
    boundary = 1
    value = 0.0
  [../]
  [./bottom_fix_y]
    type = PresetBC
    variable = disp_y
    boundary = 1
    value = 0.0
  [../]
  [./bottom_fix_z]
    type = PresetBC
    variable = disp_z
    boundary = 1
    value = 0.0
  [../]

  [./convective]
    type = ConvectiveFluxFunction # Convective flux, e.g. q'' = h*(Tw - Tf)
    variable = temp
    boundary = 'moving_boundary'
    coefficient = 1e-1
    T_infinity = 300
    marker_uo = activated_elem_uo
  [../]
  [./convective_substrate]
    type = ConvectiveFluxFunction # Convective flux, e.g. q'' = h*(Tw - Tf)
    variable = temp
    boundary = 2
    coefficient = 1e-1
    T_infinity = 300
  [../]
[]

[Materials]
  [./youngs_modulus]
    type = ADPiecewiseLinearInterpolationMaterial
    xy_data = '0          10e+8
               599.9999   10e+8
               600        9.94e+6
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
  [./strain]
    type = ADComputeFiniteStrain
    eigenstrain_names = eigenstrain
  [../]
  [./stress]
    type = ADComputeFiniteStrainElasticStress
  [../]
  [./thermal]
    type = ADComputeThermalExpansionEigenstrain
    stress_free_temperature = 300
    thermal_expansion_coeff = 2e-8
    temperature = temp
    eigenstrain_name = eigenstrain
  [../]
  [./heat]
    type = ADHeatConductionMaterial
    specific_heat = 603
    thermal_conductivity = 10e-2
  [../]
  [./volumetric_heat]
    type = DoubleEllipsoidHeatSource
    a = 5
    b = 5
    c = 3
    power = 1000
    efficienty = 0.5
    factor = 2
    velocity = -8.47
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
