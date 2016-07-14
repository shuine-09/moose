[GlobalParams]
  order = FIRST
  family = LAGRANGE
[]

[XFEM]
  #             x0   y0   x1   y1   t0    t1
  cut_data = ' 0.0  0.5  1.0  0.5  0.0   0.0'
  qrule = volfrac
  output_cut_plane = true
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 3
  ny = 3
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
[]

[SolidMechanics]
  [./solid]
    disp_x = disp_x
    disp_y = disp_y
    use_displaced_mesh = true
  [../]
[]

[Functions]
  [./pull]
    type = PiecewiseLinear
    x='0  1   2'
    y='0  0.01 0.02'
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

[Constraints]
  [./xfem_constraint_disp_x]
    type = XFEMDisplacementConstraint
    variable = disp_x
    alpha = 1e7
    use_displaced_mesh = true
 [../]
 [./xfem_constraint_disp_y]
    type = XFEMDisplacementConstraint
    variable = disp_y
    alpha = 1e7
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
#  formulation = NonlinearPlaneStrain
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
  nl_abs_tol = 1e-10

# time control
  start_time = 0.0
  dt = 1.0
  end_time = 2.0
  num_steps = 5000

  max_xfem_update = 1
[]

[Outputs]
  file_base = penalty_glued_contact_out
  exodus = true
  execute_on = timestep_end
  [./console]
    type = Console
    perf_log = true
    output_linear = true
  [../]
[]
