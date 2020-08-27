# Input file for an oxide growing on top of a zirconium nuclear fuel cladding
# using the C4 model to compute the growth rate.
# The variable is the reduced oxygen concentration [/um^3] over Czr.
# The length unit is the micrometer.
# There's 1 moving interface (alpha/oxide).
# The ICs are the simplest "const" IC (profile constant in each phase). No gradients as input for ICs.

[GlobalParams]
  order = FIRST
  family = LAGRANGE
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 21
  ny = 20
  xmin = 0
  xmax = 600
  ymin = 0
  ymax = 600
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
    cut_data = '590 0 590 600 0 0'
    heal_always = true
    interface_velocity = velocity_ox_a
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
    function = 'if(x<590,0.0075,0.4547)'
  [../]
[]

[AuxVariables]
  [./ls_ox_a]
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

  start_time = 0
  dt = 10
  num_steps = 150
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
