[GlobalParams]
  order = FIRST
  family = LAGRANGE
[]

[XFEM]
  #             x0   y0   x1   y1   t0    t1
  cut_data = ' 1.0  0.5  0.2  0.5  0.0   10.0'
  qrule = volfrac
  output_cut_plane = true
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 41
  ny = 41
  xmin = 0.0
  xmax = 1.0
  ymin = 0.0
  ymax = 1.0
  elem_type = QUAD4
  displacements = 'disp_x disp_y'
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./u]
    initial_condition = 1
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
  [./hydrostatic]
    order = FIRST
    family = MONOMIAL
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
  [./stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
  [../]
  [./stress_xy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xy
    index_i = 0
    index_j = 1
  [../]
  [./hydrostatic]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = hydrostatic
    scalar_type = hydrostatic
  [../]
[] # AuxKernels

[Kernels]
  [./TensorMechanics]
    displacements = 'disp_x disp_y'
  [../]
[]

[Kernels]
  [./diff]
    type = DiffusionStress
    variable = u
    pressure = hydrostatic
    coefficent = 1.0e-4
  [../]
  [./time_u]
    type = TimeDerivative
    variable = u
  [../]
[]

[Functions]
  [./pull]
    type = PiecewiseLinear
    x='0  50 100'
    y='0  0.02 1'
  [../]
[]

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
    type = PresetBC
    boundary = top
    variable = disp_x
    value = 0.0
  [../]
  [./topy]
    type = FunctionPresetBC
    boundary = top
    variable = disp_y
    function = pull
  [../]
[]

[Materials]
  [./Elasticity_tensor]
    type = ComputeElasticityTensor
    block = 0
    fill_method = symmetric_isotropic
    C_ijkl = '0.3 1e6'
  [../]
  [./strain]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
    block = 0
  [../]
  [./stress]
    type = ComputeLinearElasticStress
    block = 0
  [../]
[]

[Executioner]
  type = Transient

  solve_type = 'PJFNK'
#  petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type -pc_hypre_boomeramg_max_iter'
#  petsc_options_value = '201                hypre    boomeramg      8'

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu     superlu_dist'

  line_search = 'none'

  [./Predictor]
    type = SimplePredictor
    scale = 1.0
  [../]

# controls for linear iterations
  l_max_its = 100
  l_tol = 1e-2

# controls for nonlinear iterations
  nl_max_its = 15
  nl_rel_tol = 1e-9
  nl_abs_tol = 1e-9

# time control
  start_time = 0.0
  dt = 1.0
  end_time = 10.0
  num_steps = 5000

  max_xfem_update = 1
[]

[Outputs]
  exodus = true
  execute_on = timestep_end
  [./console]
    type = Console
    perf_log = true
    output_linear = true
  [../]
[]
