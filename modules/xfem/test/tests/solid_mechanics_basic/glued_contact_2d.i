[GlobalParams]
  order = FIRST
  family = LAGRANGE
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
  [./vonmises]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./elastic_strain_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./resid_y]
   # block = 0
  [../]
  [./resid_x]
  [../]
  [./traction_y]
  [../]
[]


[SolidMechanics]
  [./solid]
    disp_x = disp_x
    disp_y = disp_y
    save_in_disp_y = resid_y
    save_in_disp_x = resid_x
#    use_displaced_mesh = true
  [../]
[]

[AuxKernels]
  [./stress_xx]
    type = MaterialTensorAux
    variable = stress_xx
    tensor = stress
    index = 0
  [../]
  [./stress_yy]
    type = MaterialTensorAux
    variable = stress_yy
    tensor = stress
    index = 1
  [../]
  [./vonmises]
    type = MaterialTensorAux
    tensor = stress
    variable = vonmises
    quantity = vonmises
    execute_on = timestep_end
  [../]
  [./elastic_strain_yy]
    type = MaterialTensorAux
    variable = elastic_strain_yy
    tensor = elastic_strain
    index = 1
  [../]
[]

[Constraints]
  [./dispx_constraint]
    type = XFEMSingleVariableConstraint
    use_displaced_mesh = false
    variable = disp_x
    alpha = 10000
  [../]
  [./dispy_constraint]
    type = XFEMSingleVariableConstraint
    use_displaced_mesh = false
    variable = disp_y
    alpha = 10000
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
  [./linelast]
    type = LinearIsotropicMaterial
    block = 0
    disp_x = disp_x
    disp_y = disp_y
    poissons_ratio = 0.31
    youngs_modulus = 200.
    thermal_expansion = 0.02
    t_ref = 0.5
  [../]
[]

[Postprocessors]
  [./react_y_top]
    type = NodalSum
    variable = resid_y
    boundary = top
  [../]
  [./react_x_top]
    type = NodalSum
    variable = resid_x
    boundary = top
  [../]
  [./react_y_bottom]
    type = NodalSum
    variable = resid_y
    boundary = bottom
  [../]
  [./react_x_bottom]
    type = NodalSum
    variable = resid_x
    boundary = bottom
  [../]
  [./react_y_left]
    type = NodalSum
    variable = resid_y
    boundary = left
  [../]
  [./react_x_left]
    type = NodalSum
    variable = resid_x
    boundary = left
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
