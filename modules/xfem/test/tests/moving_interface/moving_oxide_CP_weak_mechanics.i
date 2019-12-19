# test for an oxide growing on top of a ziroconium nuclear fuel cladding
# using the Cathcart-Pawel model to compute the growth rate

[GlobalParams]
  order = FIRST
  family = LAGRANGE
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 21
  ny = 20
  xmin = 0
  xmax = 6e-4
  ymin = 0
  ymax = 1e-3
  elem_type = QUAD4
[]

[MeshModifiers]
  [./left_bottom]
    type = AddExtraNodeset
    new_boundary = 'left_bottom'
    coord = '0.0 0.0'
  [../]
[]

[XFEM]
  qrule = volfrac
  output_cut_plane = true
[]

[UserObjects]
  [./velocity]
    type = XFEMCathcartPawelOxideVelocity
    diffusivity_at_positive_level_set = 1e-5
    diffusivity_at_negative_level_set = 1e-11
    equilibrium_concentration_jump = 0.3689
    value_at_interface_uo = value_uo
  [../]
  [./value_uo]
    type = PointValueAtXFEMInterface
    variable = 'u'
    geometric_cut_userobject = 'moving_line_segments'
    execute_on = 'nonlinear'
    level_set_var = ls
  [../]
  [./moving_line_segments]
    type = MovingLineSegmentCutSetUserObject
    cut_data = '5.9e-4 0 5.9e-4 1e-3 0 0'
    heal_always = true
    interface_velocity = velocity
  [../]
[]

[Variables]
  [./u]
  [../]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
[]

[ICs]
  [./ic_u]
    type = FunctionIC
    variable = u
    function = 'if(x<5.9e-4, 0.03, 0.6667)'
  [../]
[]

[AuxVariables]
  [./temp]
    initial_condition = 1200
  [../]
  [./ls]
    order = FIRST
    family = LAGRANGE
  [../]
  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./a_strain_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./a_strain_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./a_strain_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./b_strain_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./b_strain_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./b_strain_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

# Need to use XFEMTwoSideDirichlet that has been removed
[Constraints]
  [./u_constraint]
    type = XFEMTwoSideDirichlet
    geometric_cut_userobject = 'moving_line_segments'
    use_displaced_mesh = false
    variable = u
    value_at_positive_level_set_interface = 1.0 #0.2978
    value_at_negative_level_set_interface = 1.0 #0.6667
#    value = 0.6667
    alpha = 1e5
  [../]
  [./dispx_constraint]
    type = XFEMSingleVariableConstraint
    use_displaced_mesh = false
    variable = disp_x
    alpha = 1e8
    geometric_cut_userobject = 'moving_line_segments'
  [../]
  [./dispy_constraint]
    type = XFEMSingleVariableConstraint
    use_displaced_mesh = false
    variable = disp_y
    alpha = 1e8
    geometric_cut_userobject = 'moving_line_segments'
  [../]
[]

[Kernels]
  [./diff]
    type = ConcentrationDiffusion
    variable = u
    diffusion_coefficient_name = 'diffusion_coefficient'
  [../]
  [./time]
    type = TimeDerivative
    variable = u
  [../]
  [./TensorMechanics]
    displacements = 'disp_x disp_y'
    use_displaced_mesh = false
  [../]
[]

[AuxKernels]
  [./ls]
    type = LineSegmentLevelSetAux
    line_segment_cut_set_user_object = 'moving_line_segments'
    variable = ls
  [../]
  [./stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 0
    variable = stress_xx
  [../]
  [./stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
    variable = stress_yy
  [../]
  [./stress_xy]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 1
    variable = stress_xy
  [../]
  [./a_strain_xx]
    type = RankTwoAux
    rank_two_tensor = A_total_strain
    index_i = 0
    index_j = 0
    variable = a_strain_xx
  [../]
  [./a_strain_yy]
    type = RankTwoAux
    rank_two_tensor = A_total_strain
    index_i = 1
    index_j = 1
    variable = a_strain_yy
  [../]
  [./a_strain_xy]
    type = RankTwoAux
    rank_two_tensor = A_total_strain
    index_i = 0
    index_j = 1
    variable = a_strain_xy
  [../]
  [./b_strain_xx]
    type = RankTwoAux
    rank_two_tensor = B_total_strain
    index_i = 0
    index_j = 0
    variable = b_strain_xx
  [../]
  [./b_strain_yy]
    type = RankTwoAux
    rank_two_tensor = B_total_strain
    index_i = 1
    index_j = 1
    variable = b_strain_yy
  [../]
  [./b_strain_xy]
    type = RankTwoAux
    rank_two_tensor = B_total_strain
    index_i = 0
    index_j = 1
    variable = b_strain_xy
  [../]
[]

[Materials]
  [./diffusivity_A]
    type = GenericConstantMaterial
    prop_names = A_diffusion_coefficient
    prop_values = 1e-5
  [../]
  [./diffusivity_B]
    type = GenericConstantMaterial
    prop_names = B_diffusion_coefficient
    prop_values = 1e-11
  [../]
  [./diff_combined]
    type = LevelSetBiMaterialReal
    levelset_positive_base = 'A'
    levelset_negative_base = 'B'
    level_set_var = ls
    prop_name = diffusion_coefficient
  [../]
  [./elasticity_tensor_A]
    type = ComputeIsotropicElasticityTensor
    base_name = A
    youngs_modulus = 1e3
    poissons_ratio = 0.3
  [../]
  [./strain_A]
    type = ComputeSmallStrain
    base_name = A
    displacements = 'disp_x disp_y'
    eigenstrain_names = 'eigenstrain_A'
  [../]
  [./thermal_expansion_strain_A]
    type = ComputeThermalExpansionEigenstrain
    block = 0
    stress_free_temperature = 298
    thermal_expansion_coeff = 1.0e-5
    temperature = temp
    eigenstrain_name = A_eigenstrain_A
  [../]
  [./stress_A]
    type = ComputeLinearElasticStress
    base_name = A
  [../]
  [./elasticity_tensor_B]
    type = ComputeIsotropicElasticityTensor
    base_name = B
    youngs_modulus = 1e2
    poissons_ratio = 0.3
  [../]
  [./strain_B]
    type = ComputeSmallStrain
    base_name = B
    displacements = 'disp_x disp_y'
    eigenstrain_names = 'eigenstrain_B'
  [../]
  [./thermal_expansion_strain_B]
    type = ComputeThermalExpansionEigenstrain
    block = 0
    stress_free_temperature = 298
    thermal_expansion_coeff = 1.0e-6
    temperature = temp
    eigenstrain_name = eigenstrain_B
  [../]
  [./stress_B]
    type = ComputeLinearElasticStress
    base_name = B
  [../]
  [./combined_stress]
    type = LevelSetBiMaterialRankTwo
    levelset_positive_base = 'A'
    levelset_negative_base = 'B'
    level_set_var = ls
    prop_name = stress
  [../]
  [./combined_dstressdstrain]
    type = LevelSetBiMaterialRankFour
    levelset_positive_base = 'A'
    levelset_negative_base = 'B'
    level_set_var = ls
    prop_name = Jacobian_mult
  [../]
[]

[BCs]
# Define boundary conditions
  [./left_u]
    type = DirichletBC
    variable = u
    value = 0.0450 #0.03/0.6667
    boundary = left
  [../]

  [./right_u]
    type = DirichletBC
    variable = u
    value = 2.2388 #0.6667/0.2978
    boundary = right
  [../]

  [./bottomy]
    type = PresetBC
    boundary = left_bottom
    variable = disp_y
    value = 0.0
  [../]
  [./left_x]
    type = FunctionPresetBC
    boundary = left
    variable = disp_x
    function = 0.0
  [../]
[]

[Postprocessors]
  [./cut_data_x]
    type = OxideLayerThickness
    moving_line_segments = moving_line_segments
    execute_on = TIMESTEP_END
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  line_search = 'none'

  l_tol = 1e-3
  nl_max_its = 15
  nl_rel_tol = 1e-7
  nl_abs_tol = 1e-7

  start_time = 0.0
  dt = 10
  num_steps = 150
  max_xfem_update = 1
[]


[Outputs]
  execute_on = timestep_end
  exodus = true
  [./console]
    type = Console
    output_linear = true
  [../]
  csv = true
  perf_graph = true
[]
