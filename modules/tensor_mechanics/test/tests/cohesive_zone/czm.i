[Mesh]
  # Comment
  type = CohesiveZoneMesh
  file = czm0.e
  displacements = 'disp_x disp_y'
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
    displacements = 'disp_x disp_y'
  [../]
[]

[BCs]
  [./bottom_x]
    type = DirichletBC
    variable = disp_x
    boundary = 1
    value = 0.0
  [../]
  [./bottom_y]
    type = DirichletBC
    variable = disp_y
    boundary = 1
    value = 0.0
  [../]
  [./top_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 2
    function = 0.01*t
  [../]
  [./top_x]
    type = DirichletBC
    variable = disp_x
    boundary = 2
    value = 0
  [../]
[]

[InterfaceKernels]
  [./interface]
    type = InterfaceCohesiveZone
    variable = disp_y
    neighbor_var = disp_y
    boundary = 100
  [../]
[]

[Materials]
  [./Elasticity_tensor]
    type = ComputeElasticityTensor
    block = '1 2'
    fill_method = symmetric_isotropic
    C_ijkl = '0 0.5e6'
  [../]
  [./strain]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    block = '1 2'
  [../]
  [./stress]
    type = ComputeLinearElasticStress
    block = '1 2'
  [../]
  [./gap]
    type = MaximumNormalSeparation
    boundary = 100
    disp_x = disp_x
    disp_y = disp_y
    #block = '1 2'
  [../]
[]

[Executioner]
  # Preconditioned JFNK (default)
  type = Transient

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu     superlu_dist'

  solve_type = PJFNK
  nl_abs_tol = 1e-10
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
