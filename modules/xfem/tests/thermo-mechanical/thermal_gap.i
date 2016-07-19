[GlobalParams]
  order = FIRST
  family = LAGRANGE
[]

[XFEM]
  #             x0   y0   x1   y1   t0    t1
  cut_data = ' 0.0  1.0  1.0  1.0  0.0   0.0'
  qrule = volfrac
  output_cut_plane = true
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 21
  xmin = 0.0
  xmax = 1.0
  ymin = 0.0
  ymax = 2.0
  elem_type = QUAD4
  displacements = 'disp_x disp_y'
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./temp]
    initial_condition = 200.0     # set initial temp to ambient
  [../]
[]

[SolidMechanics]
  [./solid]
    disp_x = disp_x
    disp_y = disp_y
    use_displaced_mesh = true
  [../]
[]

[Kernels]
  [./heat]         # gradient term in heat conduction equation
    type = HeatConduction
    variable = temp
    use_displaced_mesh = true
    block = 0
  [../]
  [./heat_source]
    type = XFEMHeatSource
    variable = temp
    value = 1500
    function = source
    use_displaced_mesh = true
    block = 0
  [../]

[]

[Functions]
  [./pull]
    type = PiecewiseLinear
    x='0  0.1   1.0'
    y='0.0  0.025 0.025'
  [../]
  [./source]
  type = PiecewiseLinear
  x = '0 1'
  y = '0 1'
  [../]

[]

[BCs]
  [./bottom_temp]
    type = DirichletBC
    boundary = 'left right'
    variable = temp
    value = 200
  [../]
  [./top_temp]
    type = PresetBC
    boundary = 'top'
    variable = temp
    value = 200
  [../]
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

[Constraints]
  [./xfem_constraint]
    type = XFEMHeatTransferConstraint
    variable = temp
    use_displaced_mesh = true
  [../]
[]

[Materials]
  [./elast]
    type = Elastic
    block = 0
    disp_x = disp_x
    disp_y = disp_y
    poissons_ratio = 0.3
    youngs_modulus = 1e6
    thermal_expansion = 0.0
#  formulation = NonlinearPlaneStrain
  [../]
  [./thermal]
    type = HeatConductionMaterial
    block = 0
    temp = temp
    thermal_conductivity = 1.0
  [../]
  [./heat_cond]
    type = XFEMGapFluxMaterial 
    block = '0'
    temp = temp
    disp_x = disp_x
    disp_y = disp_y
    heat_transfer_coef = 1
    dirac = true
  [../]
[]

[Executioner]
  type = Transient

  solve_type = 'PJFNK'
  petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type -pc_hypre_boomeramg_max_iter'
  petsc_options_value = '201                hypre    boomeramg      8'

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
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-9

# time control
  start_time = 0.0
  dt = 0.1
  end_time = 0.85
  num_steps = 10

  max_xfem_update = 1
[]

[Outputs]
  file_base = thermal_gap
  exodus = true
  execute_on = timestep_end
  [./console]
    type = Console
    perf_log = true
    output_linear = true
  [../]
[]
