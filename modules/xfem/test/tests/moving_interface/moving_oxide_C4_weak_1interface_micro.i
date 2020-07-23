# Test for an oxide growing on top of a zirconium nuclear fuel cladding
# using the C4 model to compute the growth rate
# The variable is the reduced concentration [/um^3] over Czr
# The length unit is the micrometer
# There's 1 moving interface (alpha/oxide)
# The ICs are enforced using steady state solving at first time step
# The ICs are enforced using 1 additional interface at the first time step (in alpha)

[GlobalParams]
  order = FIRST
  family = LAGRANGE
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 601
  ny = 40
  xmin = 0
  xmax = 600
  ymin = 0
  ymax = 40
  elem_type = QUAD4
[]

[XFEM]
  qrule = volfrac
  output_cut_plane = true
[]

[UserObjects]
  [./velocity_ox_a]
    type = XFEMC4VelocityOxideWeakMicro
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
    cut_data = '599.2 0 599.2 40 0 0'
    heal_always = true
    interface_velocity = velocity_ox_a
  [../]
  [./velocity_a0]
    type = XFEMC4VelocityMetalWeak
    diffusivity_alpha = 10
    diffusivity_beta = 10
    value_at_interface_uo = value_uo_a0
  [../]
  [./value_uo_a0]
    type = PointValueAtXFEMInterface
    variable = 'u'
    geometric_cut_userobject = 'moving_line_segments_a0'
    execute_on = 'nonlinear'
    level_set_var = ls_a0
  [../]
  [./moving_line_segments_a0]
    type = MovingLineSegmentCutSetUserObject
    cut_data = '597.7 0 597.7 40 0 0'
    heal_always = true
    interface_velocity = velocity_a0
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
    function = 'if(x<599.2,if(x<597.7,0.0075,0.0075+(x-597.7)*0.2777),0.4241)'
  [../]
[]

[AuxVariables]
  [./ls_ox_a]
    order = FIRST
    family = LAGRANGE
  [../]
  [./ls_a0]
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
    value = 0.4241
    alpha = 1e5
  [../]
  [./u_constraint_alpha]
    type = XFEMEqualValueAtInterface
    geometric_cut_userobject = 'moving_line_segments_a0'
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
  [./ls_a0]
    type = LineSegmentLevelSetAux
    line_segment_cut_set_user_object = 'moving_line_segments_a0'
    variable = ls_a0
  [../]
[]

[Materials]
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
    type = LevelSetBiMaterialReal
    levelset_negative_base = 'alpha'
    levelset_positive_base = 'oxide'
    level_set_var = ls_ox_a
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
    value = 0.4547
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

  start_time = 0.1
  dt = 1
  num_steps = 500
  max_xfem_update = 1

[]

[Controls]
  [./steady]
    type = TimePeriod
    disable_objects = 'Kernels::time '
    enable_objects = 'UserObjects::velocity_a0 UserObjects::value_uo_a0 UserObjects::moving_line_segments_a0 Constraints::u_constraint_alpha AuxKernels::ls_a0'
    start_time = '0.1'
    end_time = '1.1'
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
