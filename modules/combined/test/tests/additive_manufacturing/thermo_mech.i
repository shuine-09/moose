[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  type = MeshGeneratorMesh
  uniform_refine = 0
[]

[MeshGenerators]
  [./mesh]
    type = GeneratedMeshGenerator
    nx = 6
    ny = 6
    nz = 76
    xmin = -1.5
    xmax = 1.5
    ymin = -1.5
    ymax = 1.5
    zmin = 0
    zmax = 38
    dim = 3
  [../]
  [./add_bnd]
    type = SideSetsFromAllElementFaces
    input = mesh
    block = '0'
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
    melt_temperature = 600
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
    melt_temperature = 600
    temp_aux = temp_aux
    execute_on = timestep_begin
  [../]
[]

[BCs]
  [./bottom_fix_x]
    type = PresetBC
    variable = disp_x
    boundary = back
    value = 0.0
  [../]
  [./bottom_fix_y]
    type = PresetBC
    variable = disp_y
    boundary = back
    value = 0.0
  [../]
  [./bottom_fix_z]
    type = PresetBC
    variable = disp_z
    boundary = back
    value = 0.0
  [../]

  [./convective]
    type = ConvectiveFluxFunction # Convective flux, e.g. q'' = h*(Tw - Tf)
    variable = temp
    boundary = 'moving_boundary'
    coefficient = 100
    T_infinity = 300
    marker_uo = activated_elem_uo
  [../]
[]

[Materials]
  [./youngs_modulus]
    type = ADPiecewiseLinearInterpolationMaterial
    xy_data = '0          10e+6
               599.9999   10e+6
               600        9.94e+6
               99900      10e3'
    property = youngs_modulus
    variable = temp
  [../]

  [./elasticity_tensor]
    type = ADComputeVariableIsotropicElasticityTensor
    youngs_modulus = youngs_modulus
    poissons_ratio = 0.3
  [../]
  [./strain]
    type = ADComputeFiniteStrain
#eigenstrain_names = eigenstrain
  [../]
  [./stress]
    type = ADComputeFiniteStrainElasticStress
  [../]
  [./thermal]
    type = ADComputeThermalExpansionEigenstrain
    stress_free_temperature = 300
    thermal_expansion_coeff = 1.3e-5
    temperature = temp
    eigenstrain_name = eigenstrain
  [../]
  [./heat]
    type = ADHeatConductionMaterial
    specific_heat = 603
    thermal_conductivity = 10e-3
  [../]
  [./volumetric_heat]
    type = DoubleEllipsoidHeatSource
    a = 1.5
    b = 1.5
    c = 1.5
    power = 425
    efficienty = 0.45
    factor = 1
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

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-8
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
