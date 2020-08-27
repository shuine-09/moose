# Test for an oxide growing on top of a zirconium nuclear fuel cladding
# using the C4 model to compute the growth rate
# The variable is the reduced concentration [/um^3] over Czr
# The length unit is the micrometer
# there's 2 moving interfaces (alpha/oxide and alpha/beta)
# The ICs are set as constants in each phase through ICs, no steady state
# Temperature dependence is included. No heat equation yet. Homogeneous T.


[GlobalParams]
  order = FIRST
  family = LAGRANGE
  temperature = 1473.15
[]

[Mesh]
  [./cmg]
    type = CartesianMeshGenerator
    dim = 2
    dx = '300 300'
    dy = '60'
    ix = '100 301'
    iy = '20'
  [../]
[]

[XFEM]
  qrule = volfrac
  output_cut_plane = true
[]

[UserObjects]
  [./velocity_ox_a]
    type = XFEMC4VelocityZrOxA
    value_at_interface_uo = value_uo_ox_a
  [../]
  [./value_uo_ox_a]
    type = PointValueAtXFEMInterface
    variable = 'u'
    geometric_cut_userobject = 'moving_line_segments_ox_a'
    execute_on = 'nonlinear'
    level_set_var = ls_ox_a
  [../]
  [./moving_line_segments_ox_a]
    type = MovingLineSegmentCutSetUserObject
    cut_data = '500 0 500 60 0 0'
    is_C4 = true
    oxa_interface = true
    heal_always = true
    interface_velocity = velocity_ox_a
  [../]
  [./velocity_a_b]
    type = XFEMC4VelocityZrAB
    value_at_interface_uo = value_uo_a_b
  [../]
  [./value_uo_a_b]
    type = PointValueAtXFEMInterface
    variable = 'u'
    geometric_cut_userobject = 'moving_line_segments_a_b'
    execute_on = 'nonlinear'
    level_set_var = ls_a_b
  [../]
  [./moving_line_segments_a_b]
    type = MovingLineSegmentCutSetUserObject
    cut_data = '400 0 400 60 0 0'
    is_C4 = true
    ab_interface = true
    heal_always = true
    interface_velocity = velocity_a_b
  [../]
[]

[Variables]
  [./u]
  [../]
[]

[ICs]
  [./ic_u]
    type = C4ZrICConst
    variable = u
  [../]
[]

[AuxVariables]
  [./ls_ox_a]
    order = FIRST
    family = LAGRANGE
  [../]
  [./ls_a_b]
    order = FIRST
    family = LAGRANGE
  [../]
[]


[Constraints]
  [./u_constraint_ox_a]
    type = XFEMEqualValueAtInterfaceC4aox
    geometric_cut_userobject = 'moving_line_segments_ox_a'
    use_displaced_mesh = false
    variable = u
    alpha = 1e5
  [../]
  [./u_constraint_a_b]
    type = XFEMEqualValueAtInterfaceC4ab
    geometric_cut_userobject = 'moving_line_segments_a_b'
    use_displaced_mesh = false
    variable = u
    alpha = 1e5
  [../]
[]

[Kernels]
  [./diff]
    type = MatDiffusion
    variable = u
    diffusivity = 'diffusion_coefficient'
  [../]
  [./time]
    type = TimeDerivative
    variable = u
  [../]
[]

[AuxKernels]
  [./ls_ox_a]
    type = LineSegmentLevelSetAux
    line_segment_cut_set_user_object = 'moving_line_segments_ox_a'
    variable = ls_ox_a
  [../]
  [./ls_a_b]
    type = LineSegmentLevelSetAux
    line_segment_cut_set_user_object = 'moving_line_segments_a_b'
    variable = ls_a_b
  [../]
[]


[Materials]
  [./diffusivity_beta]
    type = C4DiffusionCoefBeta
    prop_names = beta_diffusion_coefficient
  [../]
  [./diffusivity_alpha]
    type = C4DiffusionCoefAlpha
    prop_names = alpha_diffusion_coefficient
  [../]
  [./diffusivity_oxide]
    type = GenericConstantMaterial
    prop_names = oxide_diffusion_coefficient
    prop_values = 10e6
  [../]
  [./diff_combined]
    type = LevelSetTriMaterialReal
    levelset_neg_neg_base = 'beta'
    levelset_pos_neg_base = 'alpha'
    levelset_pos_pos_base = 'oxide'
    ls_var_1 = ls_a_b
    ls_var_2 = ls_ox_a
    prop_name = diffusion_coefficient
    outputs = exodus
  [../]
[]

[BCs]
# Define boundary conditions
  [./left_u]
    type = DirichletBC
    variable = u
    value = 0.0075
    boundary = left
  [../]

  [./right_u]
    type = DirichletBCRightC4Zr
    variable = u
    boundary = right
  [../]
[]

[Postprocessors]
  [./grad_a_ox]
    type = GradValueAtXFEMInterfacePostprocessor
    value_at_interface_uo = value_uo_ox_a
    side = -1
    execute_on ='initial timestep_begin final'
  [../]
  [./grad_a_b]
    type = GradValueAtXFEMInterfacePostprocessor
    value_at_interface_uo = value_uo_a_b
    side = +1
    execute_on ='initial timestep_begin final'
  [../]
  [./grad_b_a]
    type = GradValueAtXFEMInterfacePostprocessor
    value_at_interface_uo = value_uo_a_b
    side = -1
    execute_on ='initial timestep_begin final'
  [../]
  [./position_ox_a]
    type = PositionOfXFEMInterfacePostprocessor
    value_at_interface_uo = value_uo_ox_a
    execute_on ='timestep_end final'
  [../]
  [./position_a_b]
    type = PositionOfXFEMInterfacePostprocessor
    value_at_interface_uo = value_uo_a_b
    execute_on ='timestep_end final'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  line_search = 'none'



  l_tol = 1e-3
  l_max_its = 10
  nl_max_its = 15
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6

  start_time = 19
  dt = 1
  num_steps = 1481
  max_xfem_update = 1

[]


[Outputs]
  execute_on = timestep_end
  exodus = true
  [./console]
    type = Console
    output_linear = true
  [../]
  csv = true
  perf_graph = true
[]
