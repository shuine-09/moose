[GlobalParams]
  order = FIRST
  family = LAGRANGE
  displacements = 'disp_x disp_y'
  volumetric_locking_correction = true
[]

[XFEM]
  geometric_cut_userobjects = 'level_set_cut_uo'
  qrule = volfrac
  output_cut_plane = true
[]

[UserObjects]
  [./level_set_cut_uo]
    type = LevelSetCutUserObject
    level_set_var = ls
    execute_on = 'nonlinear'
  [../]
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 31
  ny = 31
  xmin = 0.0
  xmax = 5.
  ymin = 0.0
  ymax = 5.
  elem_type = QUAD4
  displacements = 'disp_x disp_y'
[]

# [MeshModifiers]
#   [./left_bottom]
#     type = AddExtraNodeset
#     new_boundary = 'left_bottom'
#     coord = '0.0 0.0'
#   [../]
#   [./left_top]
#     type = AddExtraNodeset
#     new_boundary = 'left_top'
#     coord = '0.0 5.'
#   [../]
# []

[AuxVariables]
  [./ls]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxKernels]
  [./ls_function]
    type = FunctionAux
    variable = ls
    function = ls_func
  [../]
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./temp]
    initial_condition = 600.0
  [../]
[]

[Functions]
  [./ls_func]
    type = ParsedFunction
    value = 'sqrt((y-2.5)*(y-2.5) + (x-2.5)*(x-2.5)) - 1.5'
  [../]
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
  [./a_creep_strain_yy]
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
  [./eff_creep_strain]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  [./TensorMechanics]
    use_displaced_mesh = true
    decomposition_method = EigenSolution
    temperature = temp
  [../]
  [./heat]
    type = HeatConduction
    variable = temp
  [../]
  [./heat_ie]
    type = HeatConductionTimeDerivative
    variable = temp
  [../]
[]

[AuxKernels]
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
  [./a_creep_strain_yy]
    type = RankTwoAux
    rank_two_tensor = A_creep_strain
    variable = a_creep_strain_yy
    index_i = 1
    index_j = 1
    execute_on = timestep_end
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
  [./eff_creep_strain]
    type = MaterialRealAux
    property = effective_creep_strain
    variable = eff_creep_strain
  [../]
[]

[Constraints]
  [./dispx_constraint]
    type = XFEMSingleVariableConstraint
    use_displaced_mesh = false
    variable = disp_x
    alpha = 1e8
  [../]
  [./dispy_constraint]
    type = XFEMSingleVariableConstraint
    use_displaced_mesh = false
    variable = disp_y
    alpha = 1e8
  [../]
  [./temp]
    type = XFEMSingleVariableConstraint
    use_displaced_mesh = false
    variable = temp
    alpha = 1e6
  [../]
[]

[Functions]
  [./appl_dispy]
    type = PiecewiseLinear
    x = '0     1.0     2.0'
    y = '0.0 0.25e-4 0.50e-4'
  [../]
[]

[BCs]
  [./left_x]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]
  [./bot_y]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]
  [./top_y]
    type = FunctionPresetBC
    variable = disp_y
    boundary = top
    function = appl_dispy
  [../]
  [./temp_fix]
    type = DirichletBC
    variable = temp
    boundary = 'top right'
    value = 600.0
  [../]
[]

[Materials]
  [./elasticity_tensor_A]
    type = ComputeIsotropicElasticityTensor
    base_name = A
    youngs_modulus = 250e9
    poissons_ratio = 0.25
  [../]
  # [./strain_A]
  #   type = ComputePlaneFiniteStrain
  #   base_name = A
  # [../]
  # [./stress_A]
  #   type = ComputeLinearElasticStress
  #   base_name = A
  # [../]

  [./strain_A]
    type = ComputeFiniteStrain #ComputePlaneFiniteStrain 
    block = 0
    base_name = A
  [../]

  # [./stress_A]
  #   type = ComputeFiniteStrainElasticStress
  #   base_name = A
  # [../]

  [./radial_return_stress_A]
    type = ComputeMultipleInelasticStress
    block = 0
    inelastic_models = 'powerlawcrp_A'
    base_name = A
  [../]
  [./powerlawcrp_A]
    type = PowerLawCreepStressUpdate
    block = 0
    coefficient = 3.125e-14
    n_exponent = 2.0
    m_exponent = 0.0
    activation_energy = 0.0
    creep_prepend = A_
  [../]

  [./elasticity_tensor_B]
    type = ComputeIsotropicElasticityTensor
    base_name = B
    youngs_modulus = 250e9
    poissons_ratio = 0.25
  [../]

  [./strain_B]
    type = ComputeFiniteStrain
    base_name = B
  [../]
  [./stress_B]
    type = ComputeFiniteStrainElasticStress
    base_name = B
  [../]

  [./combined]
    type = LevelSetMultiStressMaterial
    levelset_plus_base = 'B'
    levelset_minus_base = 'A'
    level_set_var = ls
  [../]

  [./thermal]
    type = HeatConductionMaterial
    block = 0
    specific_heat = 1.0
    thermal_conductivity = 100.
  [../]
  [./density]
    type = Density
    block = 0
    density = 1.0
  [../]
[]

[Executioner]
  type = Transient

  solve_type = 'PJFNK'
  #petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type -pc_hypre_boomeramg_max_iter'
  #petsc_options_value = '201                hypre    boomeramg      8'


   petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
   petsc_options_value = 'lu     superlu_dist'

   #line_search = 'bt'

  #[./Predictor]
  #  type = SimplePredictor
  #  scale = 1.0
  #[../]

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6
  l_tol = 1e-6
  l_max_its = 20
  nl_max_its = 30

  dt = 0.1
  start_time = 0.0
  num_steps = 100
  #end_time = 2.0

  max_xfem_update = 1
[]

[Outputs]
  exodus = true
  execute_on = timestep_end
  csv = true
  [./console]
    type = Console
    perf_log = true
    output_linear = true
  [../]
[]
