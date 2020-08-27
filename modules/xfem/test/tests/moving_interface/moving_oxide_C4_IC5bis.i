# Input file for an oxide growing on top of a zirconium nuclear fuel cladding
# using the C4 model to compute the growth rate.
# The variable is the reduced concentration [/um^3] over Czr.
# The length unit is the micrometer.
# There are 2 moving interfaces (alpha/oxide and alpha/beta).
# The ICs are enforced using steady state solving at first time step.
# The ICs are enforced using 3 additional interfaces at the first time step (2 in alpha, 1 in beta)
# The difference with IC5 is that the profile is constant between the two additional alpha interfaces

# Fixed T=1200C.
# File used to study influence of ICs.

[GlobalParams]
  order = FIRST
  family = LAGRANGE
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 601
  ny = 20
  xmin = 0
  xmax = 600
  ymin = 0
  ymax = 20
  elem_type = QUAD4
[]

[XFEM]
  qrule = volfrac
  output_cut_plane = true
[]

[UserObjects]
  [./velocity_ox_a]
    type = XFEMC4VelocityZrOxA
    diffusivity_alpha = 10
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
    cut_data = '590 0 590 20 0 0'
    heal_always = true
    interface_velocity = velocity_ox_a
  [../]
  [./velocity_a_b]
    type = XFEMC4VelocityZrAB
    diffusivity_alpha = 10
    diffusivity_beta = 60
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
    cut_data = '572.2 0 572.2 20 0 0'
    heal_always = true
    interface_velocity = velocity_a_b
  [../]
  [./velocity_a1]
    type = XFEMC4VelocityZrAB
    diffusivity_alpha = 10
    diffusivity_beta = 10
    value_at_interface_uo = value_uo_a1
  [../]
  [./value_uo_a1]
    type = PointValueAtXFEMInterface
    variable = 'u'
    geometric_cut_userobject = 'moving_line_segments_a1'
    execute_on = 'nonlinear'
    level_set_var = ls_a1
  [../]
  [./moving_line_segments_a1]
    type = MovingLineSegmentCutSetUserObject
    cut_data = '574.0 0 574.0 20 0 0'
    heal_always = true
    interface_velocity = velocity_a1
  [../]
  [./velocity_a2]
    type = XFEMC4VelocityZrAB
    diffusivity_alpha = 10
    diffusivity_beta = 10
    value_at_interface_uo = value_uo_a2
  [../]
  [./value_uo_a2]
    type = PointValueAtXFEMInterface
    variable = 'u'
    geometric_cut_userobject = 'moving_line_segments_a2'
    execute_on = 'nonlinear'
    level_set_var = ls_a2
  [../]
  [./moving_line_segments_a2]
    type = MovingLineSegmentCutSetUserObject
    cut_data = '579.15 0 579.15 20 0 0'
    heal_always = true
    interface_velocity = velocity_a2
  [../]
  [./velocity_b0]
    type = XFEMC4VelocityZrAB
    diffusivity_alpha = 60
    diffusivity_beta = 60
    value_at_interface_uo = value_uo_b0
  [../]
  [./value_uo_b0]
    type = PointValueAtXFEMInterface
    variable = 'u'
    geometric_cut_userobject = 'moving_line_segments_b0'
    execute_on = 'nonlinear'
    level_set_var = ls_b0
  [../]
  [./moving_line_segments_b0]
    type = MovingLineSegmentCutSetUserObject
    cut_data = '531.0 0 531.0 20 0 0'
    heal_always = true
    interface_velocity = velocity_b0
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
    function = 'if(x<590.0,0.0075,0.3679)'
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
  [./ls_a1]
    order = FIRST
    family = LAGRANGE
  [../]
  [./ls_a2]
    order = FIRST
    family = LAGRANGE
  [../]
  [./ls_b0]
    order = FIRST
    family = LAGRANGE
  [../]
[]


[Constraints]
  [./u_constraint_ox_a]
    type = XFEMEqualValueAtInterface
    geometric_cut_userobject = 'moving_line_segments_ox_a'
    use_displaced_mesh = false
    variable = u
    value = 0.3373
    alpha = 1e5
  [../]
  [./u_constraint_a_b]
    type = XFEMEqualValueAtInterface
    geometric_cut_userobject = 'moving_line_segments_a_b'
    use_displaced_mesh = false
    variable = u
    value = 0.0373
    alpha = 1e5
  [../]
  [./u_constraint_alpha1]
    type = XFEMEqualValueAtInterface
    geometric_cut_userobject = 'moving_line_segments_a1'
    use_displaced_mesh = false
    variable = u
    value = 0.0485
    alpha = 1e5
  [../]
  [./u_constraint_alpha2]
    type = XFEMEqualValueAtInterface
    geometric_cut_userobject = 'moving_line_segments_a2'
    use_displaced_mesh = false
    variable = u
    value = 0.0485
    alpha = 1e5
  [../]
  [./u_constraint_beta]
    type = XFEMEqualValueAtInterface
    geometric_cut_userobject = 'moving_line_segments_b0'
    use_displaced_mesh = false
    variable = u
    value = 0.0075
    alpha = 1e5
  [../]

[]

[Kernels]
  [./diff]
    type = ConcentrationDiffusion
    variable = u
    diffusion_coefficient_name = 'diffusion_coefficient'
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
  [./ls_a1]
    type = LineSegmentLevelSetAux
    line_segment_cut_set_user_object = 'moving_line_segments_a1'
    variable = ls_a1
  [../]
  [./ls_a2]
    type = LineSegmentLevelSetAux
    line_segment_cut_set_user_object = 'moving_line_segments_a2'
    variable = ls_a2
  [../]
  [./ls_b0]
    type = LineSegmentLevelSetAux
    line_segment_cut_set_user_object = 'moving_line_segments_b0'
    variable = ls_b0
  [../]
[]

[Materials]
  [./diffusivity_beta]
    type = GenericConstantMaterial
    prop_names = beta_diffusion_coefficient
    prop_values = 60
  [../]
  [./diffusivity_alpha]
    type = GenericConstantMaterial
    prop_names = alpha_diffusion_coefficient
    prop_values = 10
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
    type = DirichletBC
    variable = u
    value = 0.3679
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

  start_time = 20
  dt = 1
  num_steps = 480
  max_xfem_update = 1

[]

[Controls]
  [./steady]
    type = TimePeriod
    disable_objects = 'Kernels::time '
    enable_objects = 'UserObjects::velocity_a1 UserObjects::value_uo_a1 UserObjects::moving_line_segments_a1 Constraints::u_constraint_alpha1 AuxKernels::ls_a1
                      UserObjects::velocity_a2 UserObjects::value_uo_a2 UserObjects::moving_line_segments_a2 Constraints::u_constraint_alpha2 AuxKernels::ls_a2
                      UserObjects::velocity_b0 UserObjects::value_uo_b0 UserObjects::moving_line_segments_b0 Constraints::u_constraint_beta AuxKernels::ls_b0 '
    start_time = '20'
    end_time = '21'
  [../]
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
