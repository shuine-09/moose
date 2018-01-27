[GlobalParams]
  order = FIRST
  family = LAGRANGE
[]

[XFEM]
  geometric_cut_userobjects = 'line_seg_cut_uo'
  qrule = moment_fitting
  output_cut_plane = true
[]

[UserObjects]
  [./line_seg_cut_uo]
    type = LineSegmentCutUserObject
    cut_data = '0.0 5 5 5'
    #cut_data = '0.0 5. 5 5.'
    time_start_cut = 0.0
    time_end_cut = 0.0
  [../]
  [./pair_qps]
    type = XFEMElementPairQPProvider
  [../]
  [./manager]
    type = XFEMElemPairMaterialManager
    material_names = 'elasticity_tensor_interface strain_interface material1'
    element_pair_qps = pair_qps
    use_neighbor = false
  [../]
  [./manager_neighbor]
    type = XFEMElemPairMaterialManager
    material_names = 'neighbor_elasticity_tensor_interface neighbor_strain_interface neighbor_material1'
    element_pair_qps = pair_qps
    use_neighbor = true
  [../]
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 2
  ny = 3
  xmin = 0.0
  xmax = 5.
  ymin = 0.0
  ymax = 10
  elem_type = QUAD
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

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
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
[]


[Kernels]
  [./TensorMechanics]
    displacements = 'disp_x disp_y'
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
[]

[Constraints]
  [./dispx_constraint]
    type = XFEMNitscheDisplacementConstraint
    use_displaced_mesh = false
    component = 0
    variable = disp_x
    disp_x = disp_x
    disp_y = disp_y
    stiffness = 200
    max_traction = 1
    Gc = 0.005
    manager = manager
    manager_neighbor = manager_neighbor
    base_name = interface
  [../]
  [./dispy_constraint]
    type = XFEMNitscheDisplacementConstraint
    use_displaced_mesh = false
    component = 1
    variable = disp_y
    disp_x = disp_x
    disp_y = disp_y
    stiffness = 200
    max_traction = 1
    Gc = 0.005
    manager = manager
    manager_neighbor = manager_neighbor
    base_name = interface
  [../]
[]

[Functions]
  [./pull_up_and_down]
    type = ParsedFunction
    value = 'if(t < 8, 0.001 * t, if(t <= 16, -0.001 * (t-8) + 0.008, (t-16) * 0.001))'
  [../]
[]

[BCs]
  # [./bottomx]
  #   type = PresetBC
  #   boundary = left_bottom
  #   variable = disp_x
  #   value = 0.0
  # [../]
  # [./bottomy]
  #   type = PresetBC
  #   boundary = bottom
  #   variable = disp_y
  #   value = 0.0
  # [../]
  # [./topx]
  #   type = PresetBC
  #   boundary = left_top
  #   variable = disp_x
  #   value = 0.0
  # [../]
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
   type = PresetBC
   boundary = top
   variable = disp_x
   value = 0.0
  [../]
  [./topy]
#    type = PresetBC
    type = FunctionPresetBC
    boundary = top
    variable = disp_y
    function = '-0.001*t'
    #function = pull_up_and_down
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 50000
    poissons_ratio = 0.3
  [../]
  [./strain]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
  [../]
  [./stress]
    type = ComputeLinearElasticStress
  [../]
  [./material1]
    type = MaximumNormalSeparation
    compute = false
    disp_x = disp_x
    disp_y = disp_y
  [../]
  [./neighbor_material1]
    type = MaximumNormalSeparation
    compute = false
    disp_x = disp_x
    disp_y = disp_y
    base_name = neighbor
  [../]
  # [./stress_interface]
  #   type = ComputeLinearElasticStress
  #   compute = false
  #   base_name = interface
  # [../]
  [./elasticity_tensor_interface]
    type = ComputeIsotropicElasticityTensor
    compute = false
    youngs_modulus = 50000
    poissons_ratio = 0.3
    base_name = interface
  [../]
  [./neighbor_elasticity_tensor_interface]
    type = ComputeIsotropicElasticityTensor
    compute = false
    youngs_modulus = 100000
    poissons_ratio = 0.3
    base_name = neighbor_interface
  [../]
  [./strain_interface]
    type = ComputeSmallStrain
    compute = false
    displacements = 'disp_x disp_y'
    base_name = interface
  [../]
  [./neighbor_strain_interface]
    type = ComputeSmallStrain
    compute = false
    displacements = 'disp_x disp_y'
    base_name = neighbor_interface
  [../]
[]

#[UserObjects]
#  [./calc_trac_x]
#    type = ComputeCohesiveTraction
#    execute_on = timestep_end
#  [../]
#  [./calc_trac_y]
#    type = ComputeCohesiveTraction
#    execute_on = timestep_end
#  [../]
#[]

[Postprocessors]
  [./disp_y]
    type = NodalMaxValue
    variable = disp_y
    boundary = top
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

  solve_type = 'PJFNK'
  #petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type -pc_hypre_boomeramg_max_iter'
  #petsc_options_value = '201                hypre    boomeramg      8'


   petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
   petsc_options_value = 'lu     superlu_dist'

   line_search = bt

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
  dt = 1
  end_time = 40.0
  num_steps = 40

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
