[GlobalParams]
  order = FIRST
  family = LAGRANGE
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 41
  ny = 21
  xmin = 0
  xmax = 2
  ymin = 0
  ymax = 1
  elem_type = QUAD4
[]

[XFEM]
  qrule = volfrac
  output_cut_plane = true
[]

[MeshModifiers]
  [./add]
    type = BoundingBoxNodeSet
    new_boundary = 'all'
    top_right = '2 1 0'
    bottom_left = '0 0 0'
  [../]
[]

[UserObjects]
  [./moving_line_segments]
    type = MovingLineSegmentCutSetUserObject
    cut_data = '0.5 0 0.5 1.0 0 0'
    heal_always = true
    var = u
    interface_value_uo = value_uo
    diffusivity_at_positive_level_set_side = 1
    diffusivity_at_negative_level_set_side = 1
  [../]
[]

[Variables]
  [./u]
  [../]
[]

[ICs]
  [./ic_u]
    type = FunctionIC
    variable = u
    function = 'if(x<0.5, 2, 1)'
  [../]
[]

[AuxVariables]
  [./ls]
    order = FIRST
    family = LAGRANGE
  [../]
  [./u_scale]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Constraints]
  [./u_constraint]
    type = XFEMTwoSideDirichletNitsche
    geometric_cut_userobject = 'moving_line_segments'
    use_displaced_mesh = false
    variable = u
    value_at_positive_level_set_interface = 2 #10
    value_at_negative_level_set_interface = 2
    alpha = 1e5
    diffusivity_at_positive_level_set_side = 5 #1
    diffusivity_at_negative_level_set_side = 1
  [../]
[]

[Kernels]
  [./diff]
    type = ConcentrationDiffusion
    variable = u
  [../]
  [./time]
    type = TimeDerivative
    variable = u
  [../]
[]

[AuxKernels]
  [./ls]
    type = LineSegmentLevelSetAux
    geometric_cut_userobject = 'moving_line_segments'
    variable = ls
  [../]
  [./u_scale]
    type = VariableScaleAux
    scale_var = u
    variable = u_scale
    scale_factor = 5
  [../]
[]

[Materials]
  [./diffusivity_A]
    type = GenericConstantMaterial
    prop_names = A_diffusion_coefficient
    prop_values = 5
  [../]

  [./diffusivity_B]
    type = GenericConstantMaterial
    prop_names = B_diffusion_coefficient
    prop_values = 1 #1.0e-5
  [../]

  [./diff_combined]
    type = LevelSetBiMaterialProperty
    levelset_positive_base = 'A'
    levelset_negative_base = 'B'
    level_set_var = ls
    prop_name = diffusion_coefficient
  [../]

  [./time_A]
    type = GenericConstantMaterial
    prop_names = A_time_step_scale
    prop_values = 5
  [../]

  [./time_B]
    type = GenericConstantMaterial
    prop_names = B_time_step_scale
    prop_values = 1
  [../]

  [./time_combined]
    type = LevelSetBiMaterialProperty
    levelset_positive_base = 'A'
    levelset_negative_base = 'B'
    level_set_var = ls
    prop_name = time_step_scale
  [../]
[]

[BCs]
# Define boundary conditions
  [./left_u]
    type = DirichletBC
    variable = u
    value = 2
    boundary = 3
  [../]

  [./right_u]
    type = NeumannBC
    variable = u
    boundary = 1
    value = 0
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  # petsc_options_iname = '-pc_type -pc_hypre_type'
  # petsc_options_value = 'hypre boomeramg'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  line_search = 'none'

  l_tol = 1e-3
  nl_max_its = 15
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-9

  start_time = 0.0
  dt = 0.01
  end_time = 1.5
  num_steps = 10
  max_xfem_update = 1
[]

[UserObjects]
  [./value_uo]
    type = PointValueAtXFEMInterface
    variable = 'u'
    geometric_cut_userobject = 'moving_line_segments'
    execute_on = 'nonlinear'
    level_set_var = ls
  [../]
[]

[Outputs]
  execute_on = timestep_end
  exodus = true
  [./console]
    type = Console
    perf_log = true
    output_linear = true
  [../]
  csv = true
[]
