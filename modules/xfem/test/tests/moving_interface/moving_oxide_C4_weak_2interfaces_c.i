# test for an oxide growing on top of a zirconium nuclear fuel cladding
# using the C4 model to compute the growth rate
# the variable is the reduced concentration [/m^3] over Czr

[GlobalParams]
  order = FIRST
  family = LAGRANGE
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 61
  ny = 20
  xmin = 0
  xmax = 6e-4
  ymin = 0
  ymax = 1e-3
  elem_type = QUAD4
[]

[XFEM]
  qrule = volfrac
  output_cut_plane = true
[]

[UserObjects]
  [./velocity_ox_a]
    type = XFEMC4VelocityOxideWeak
    diffusivity_alpha = 1e-11
    value_at_interface_uo = value_uo_ox_a
    x0 = 5.9e-4
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
    cut_data = '5.9e-4 0 5.9e-4 1e-3 0 0'
    heal_always = true
    interface_velocity = velocity_ox_a
  [../]
  [./velocity_a_b]
    type = XFEMC4VelocityMetalWeak
    diffusivity_alpha = 1e-11
    diffusivity_beta = 6e-11
    value_at_interface_uo = value_uo_a_b
    x0 = 5.7e-4
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
    cut_data = '5.7e-4 0 5.7e-4 1e-3 0 0'
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
    type = FunctionIC
    variable = u
    function = 'if(x<5.9e-4,if(x<5.7e-4,if(x<5e-4,0.0075,0.0075+(x-5e-4)*425.71),0.0373+(x-5.7e-4)*1.5e4), 0.3679)'
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
[]

[Materials]
  [./diffusivity_beta]
    type = GenericConstantMaterial
    prop_names = beta_diffusion_coefficient
    prop_values = 6e-11
  [../]
  [./diffusivity_alpha]
    type = GenericConstantMaterial
    prop_names = alpha_diffusion_coefficient
    prop_values = 1e-11
  [../]
  [./diffusivity_oxide]
    type = GenericConstantMaterial
    prop_names = oxide_diffusion_coefficient
    prop_values = 1e-5
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

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  line_search = 'none'



  l_tol = 1e-3
  nl_max_its = 15
  nl_rel_tol = 1e-7
  nl_abs_tol = 1e-7

  start_time = 0.0
  dt = 10
  num_steps = 150.0
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
