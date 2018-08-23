[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  file = n10-id1.msh
# file = grains_debug.e
#  file = poly_3d.e
#  file = gmsh_test_in.e
#  file = hex_grains_v2.e
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
                1  0  0
                0  0  -1
                0  0  1'
     fixed_normal = true
     new_boundary = '101 102 103 104 105 106'
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
  [./disp_z]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./TensorMechanics]
  [../]
[]

[BCs]
  [./x1]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 103
    function = -0.01*t
  [../]
  [./x2]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 104
    function = 0.01*t
  [../]
  # [./bottom_y]
  #   type = DirichletBC
  #   variable = disp_y
  #   boundary = 101
  #   value = 0.0
  # [../]
  # [./bottom_z]
  #   type = DirichletBC
  #   variable = disp_y
  #   boundary = 101
  #   value = 0.0
  # [../]
  [./y1]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 101
    function = -0.01*t
  [../]
  [./y2]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 102
    function = 0.01*t
  [../]
  [./z1]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = 105
    function = -0.01*t
  [../]
  [./z2]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = 106
    function = 0.01*t
  [../]
[]

[InterfaceKernels]
  [./interface]
    type = InterfaceCohesiveZone
    variable = disp_y
    neighbor_var = disp_y
    boundary = 'interface'
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  #  ksp_norm = default
  [../]
[]

[Materials]
  [./Elasticity_tensor]
    type = ComputeElasticityTensor
    #block = '1 2'
    fill_method = symmetric_isotropic_E_nu
    C_ijkl = '0.5e6 0'
  [../]
  [./strain]
    type = ComputeSmallStrain
    #block = '1 2'
  [../]
  [./stress]
    type = ComputeLinearElasticStress
    #block = '1 2'
  [../]
  # [./gap]
  #   type = MaximumNormalSeparation
  #   boundary = 'Block1_Block2'
  #   disp_x = disp_x
  #   disp_y = disp_y
  # [../]
[]

[Executioner]
  # Preconditioned JFNK (default)
  type = Transient

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu     superlu_dist'

  solve_type = PJFNK
  nl_abs_tol = 1e-6
  l_tol = 1.0e-4
  l_max_its = 5
  start_time = 0.0
  dt = 1.0
  num_steps = 2
  end_time = 2.0
[]

[Outputs]
  [./out]
    type = Exodus
  [../]
[]
