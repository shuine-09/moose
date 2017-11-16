[GlobalParams]
  order = FIRST
  family = LAGRANGE
  displacements = 'disp_x disp_y'
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
  nx = 5
  ny = 11
  xmin = 0.0
  xmax = 5.
  ymin = 0.0
  ymax = 10
  elem_type = QUAD4
  displacements = 'disp_x disp_y'
[]

[MeshModifiers]
  [./left_bottom]
    type = AddExtraNodeset
    new_boundary = 'left_bottom'
    coord = '0.0 0.0'
  [../]
  [./left_top]
    type = AddExtraNodeset
    new_boundary = 'left_top'
    coord = '0.0 10.'
  [../]
[]

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
[]

[Functions]
  [./ls_func]
    type = ParsedFunction
    value = 'y-2.5'
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
[]

[Kernels]
  [./TensorMechanics]
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
[]

#[Constraints]
#  [./dispx_constraint]
#    type = XFEMSingleVariableConstraint
#    use_displaced_mesh = false
#    variable = disp_x
#    alpha = 10000
#  [../]
#  [./dispy_constraint]
#    type = XFEMSingleVariableConstraint
#    use_displaced_mesh = false
#    variable = disp_y
#    alpha = 10000
#  [../]
#[]

[BCs]
  [./bottomx]
    type = PresetBC
    boundary = bottom
    variable = disp_x
    value = 0.0
  [../]
  [./bottomy]
    type = PresetBC
    boundary = bottom
    variable = disp_y
    value = 0.0
  [../]
  [./topx]
    type = FunctionPresetBC
    boundary = top
    variable = disp_x
    function = 0.03*t
  [../]
  [./topy]
#    type = PresetBC
    type = FunctionPresetBC
    boundary = top
    variable = disp_y
#    value = 0.1
    function = '0.03*t'
  [../]
[]

[Materials]
  [./elasticity_tensor_A]
    type = ComputeElasticityTensor
    base_name = A
    fill_method = symmetric9
    C_ijkl = '1e6 1e5 1e5 1e6 0 1e6 .4e6 .2e6 .5e6'
  [../]
  [./strain_A]
    type = ComputeSmallStrain
    base_name = A
    eigenstrain_names = eigenstrain
  [../]
  [./stress_A]
    type = ComputeLinearElasticStress
    base_name = A
  [../]
  [./eigenstrain_A]
    type = ComputeEigenstrain
    base_name = A
    eigen_base = '0.1 0.05 0 0 0 0.01'
    prefactor = -1
    eigenstrain_name = eigenstrain
  [../]

  [./elasticity_tensor_B]
    type = ComputeElasticityTensor
    base_name = B
    fill_method = symmetric9
    C_ijkl = '1e6 0 0 1e6 0 1e6 .5e6 .5e6 .5e6'
  [../]
  [./strain_B]
    type = ComputeSmallStrain
    base_name = B
    eigenstrain_names = 'B_eigenstrain'
  [../]
  [./stress_B]
    type = ComputeLinearElasticStress
    base_name = B
  [../]
  [./eigenstrain_B]
    type = ComputeEigenstrain
    base_name = B
    eigen_base = '0.1 0.05 0 0 0 0.01'
    prefactor = -1
    eigenstrain_name = 'B_eigenstrain'
  [../]

  [./combined]
    type = LevelSetMultiStressMaterial
    phase_base = 'A  B'
    level_set_var = ls
  [../]
[]

[Executioner]
  type = Transient

  solve_type = 'PJFNK'
#  petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type -pc_hypre_boomeramg_max_iter'
#  petsc_options_value = '201                hypre    boomeramg      8'


   petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
   petsc_options_value = 'lu     superlu_dist'

   line_search = 'bt'

  #[./Predictor]
  #  type = SimplePredictor
  #  scale = 1.0
  #[../]

# controls for linear iterations
  l_max_its = 20
  l_tol = 1e-3

# controls for nonlinear iterations
  nl_max_its = 15
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-5

# time control
  start_time = 0.0
  dt = 0.1
  end_time = 1.0
  num_steps = 5000

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
