[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Mesh]
#  file = n10-id1.msh
# file = grains_debug.e
#  file = poly_3d.e
#  file = gmsh_test_in.e
  file = poly2d.msh
  parallel_type = REPLICATED
[]

[MeshModifiers]
  [./breakmesh]
    type = BreakMeshByBlock
    split_interface = false
  [../]
  [./add_side_sets]
     type = SideSetsFromNormals
     normals = '0  -1  0
                0  1  0
                -1 0  0
                1  0  0'
     fixed_normal = true
     new_boundary = '101 102 103 104'
   [../]
[]

[Variables]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./TensorMechanics]
    save_in = 'resid_x resid_y'
  [../]
[]

[Functions]
  [./pull_up_and_down]
    type = ParsedFunction
    value = 'if(t < 100, 0.0001 * t, -0.0001 * (t-100) + 0.01)'
    #value = '0.0001 * t'
  [../]
[]

[BCs]
  [./bottom_x]
    type = DirichletBC
    variable = disp_x
    boundary = 101
    value = 0.0
  [../]
  [./bottom_y]
    type = DirichletBC
    variable = disp_y
    boundary = 101
    value = 0.0
  [../]
  [./top_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 102
    function = pull_up_and_down
  [../]
  [./top_x]
    type = DirichletBC
    variable = disp_x
    boundary = 102
    value = 0
  [../]
[]

[AuxVariables]
  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./resid_y]
  [../]
  [./resid_x]
  [../]
  [./traction_y]
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
  [./stress_xy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xy
    index_i = 0
    index_j = 1
  [../]
  [./stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
  [../]
[]

[InterfaceKernels]
  [./interface_x]
    type = InterfaceCohesiveZone
    variable = disp_x
    neighbor_var = disp_x
    component = 0
    boundary = 'interface'
    stiffness = 100
    max_traction = 0.5
    Gc = 50
  [../]
  [./interface_y]
    type = InterfaceCohesiveZone
    variable = disp_y
    neighbor_var = disp_y
    component = 1
    boundary = 'interface'
    # stiffness = 10
    # max_traction = 0.05
    # Gc = 0.0005
    stiffness = 100
    max_traction = 0.5
    Gc = 50
  [../]

  [./interface_x2]
    type = InterfaceGluedContact
    variable = disp_x
    neighbor_var = disp_x
    component = 0
    boundary = 'interface'
    alpha = 10000
    tangential = true
  [../]
  [./interface_2]
    type = InterfaceGluedContact
    variable = disp_y
    neighbor_var = disp_y
    component = 1
    boundary = 'interface'
    alpha = 10000
    tangential = true
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    ksp_norm = default
  [../]
[]

[Materials]
  [./Elasticity_tensor]
    type = ComputeElasticityTensor
  #  block = '1 2'
    fill_method = symmetric_isotropic_E_nu
    C_ijkl = '5000 0.02'
  [../]
  [./strain]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
  #  block = '1 2'
  [../]
  [./stress]
    type = ComputeLinearElasticStress
  #  block = '1 2'
  [../]
  [./gap]
    type = MaximumNormalSeparation
    boundary = 'interface'
    disp_x = disp_x
    disp_y = disp_y
  [../]
[]

[Postprocessors]
  [./react_y_top]
    type = NodalSum
    variable = resid_y
    boundary = 102
  [../]
  [./react_x_top]
    type = NodalSum
    variable = resid_x
    boundary = 102
  [../]
  [./react_y_bottom]
    type = NodalSum
    variable = resid_y
    boundary = 101
  [../]
  [./react_x_bottom]
    type = NodalSum
    variable = resid_x
    boundary = 101
  [../]
  [./disp_y]
    type = NodalMaxValue
    variable = disp_y
    boundary = 102
  [../]
[]

[Executioner]
  # Preconditioned JFNK (default)
  type = Transient

  # petsc_options_iname = '-ksp_type -pc_type -pc_factor_mat_solver_package'
  # petsc_options_value = 'preonly lu     superlu_dist'

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu     superlu_dist'

  line_search = 'none'

  solve_type = PJFNK
  nl_abs_tol = 1e-8
  nl_rel_tol = 1e-8
  l_tol = 1.0e-6
  l_max_its = 5
  nl_max_its = 100
  start_time = 0.0
#  dtmin = 1
  dt = 1
  end_time = 100
  num_steps = 5000
[]

[Outputs]
  [./out]
    type = Exodus
  [../]
  csv = true
[]
