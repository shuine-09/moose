# Input file for an oxide growing outward on top of a steel 21-2N sample
# using the C4 model to compute the growth rate
# The variable is the Mn concentration [/nm^3]
# The length unit is the micrometer. The time unit is the hour
# there's 1 fixed interface (oxide/metal) and 1 moving interface (oxide/gas)
# The gas is the left part of the mesh (void)
# The ICs are set as constants in each phase through ICs, no steady state
# Homogeneous T=700C for now.

[GlobalParams]
  order = FIRST
  family = LAGRANGE
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 541
  ny = 100
  xmin = 0
  xmax = 108
  ymin = 0
  ymax = 20
  elem_type = QUAD4
[]

[XFEM]
  qrule = volfrac
  output_cut_plane = true
[]

[UserObjects]
  [./fixed_cut_oxide_metal]
    type = LineSegmentCutSetUserObject
    cut_data = '105 0 105 20 0 0'
  [../]
  [./value_uo_m_ox]
    type = PointValueAtXFEMInterface
    variable = 'C_Mn'
    geometric_cut_userobject = 'fixed_cut_oxide_metal'
    execute_on = 'nonlinear'
    level_set_var = ls_metal_oxide
  [../]
  [./moving_cut_oxide]
    type = MovingLineSegmentCutSetUserObject
    cut_data = '105.3 0 105.3 20 0 0'
    heal_always = true
    interface_velocity = velocity_oxide
  [../]
  [./velocity_oxide]
    type = XFEMVelocitySteelOxLimDiff
    value_at_interface_uo = value_uo_m_ox
  [../]
  [./value_uo_oxide]
    type = PointValueAtXFEMInterface
    variable = 'C_Mn'
    geometric_cut_userobject = 'moving_cut_oxide'
    execute_on = 'nonlinear'
    level_set_var = ls_oxide_gas
  [../]
  [./fixed_cut_grad]
    type = LineSegmentCutSetUserObject
    cut_data = '102 0 102 20 0 0'
  [../]
  [./value_uo_grad]
    type = PointValueAtXFEMInterface
    variable = 'C_Mn'
    geometric_cut_userobject = 'fixed_cut_grad'
    execute_on = 'nonlinear'
    level_set_var = ls_grad
  [../]
[]

[Variables]
  [./C_Mn]
  [../]
[]

[ICs]
  [./ic_Mn]
    type = FunctionIC
    variable = C_Mn
    function = 'if(x<105,7.1445,if(x<105.3,13.3179,13.3152))' #if(x<102,7.1445,7.1445-7.1445/3*(x-102))
  [../]
[]

[AuxVariables]
  [./ls_oxide_gas]
    order = FIRST
    family = LAGRANGE
  [../]
  [./ls_metal_oxide]
    order = FIRST
    family = LAGRANGE
  [../]
  [./ls_grad]
    order = FIRST
    family = LAGRANGE
  [../]
[]


[Constraints]
#  [./oxide_metal_constraint]
#    type = XFEMTwoSideDirichlet
#    geometric_cut_userobject = 'fixed_cut_oxide_metal'
#    use_displaced_mesh = false
#    variable = C_Mn
#    value_at_positive_level_set_interface = 0
#    value_at_negative_level_set_interface = 13.3179
#    alpha = 1e5
#  [../]
  [./oxide_metal_constraint_init]
    type = XFEMEqualValueAtInterface
    geometric_cut_userobject = 'fixed_cut_oxide_metal'
    use_displaced_mesh = false
    variable = C_Mn
    value = 0
    alpha = 1e5
  [../]
  [./oxide_gas_constraint]
    type = XFEMEqualValueAtInterface
    geometric_cut_userobject = 'moving_cut_oxide'
    use_displaced_mesh = false
    variable = C_Mn
    value = 13.3152
    alpha = 1e5
  [../]
  [./grad_constraint]
    type = XFEMEqualValueAtInterface
    geometric_cut_userobject = 'fixed_cut_grad'
   use_displaced_mesh = false
    variable = C_Mn
    value = 7.1445
    alpha = 1e5
  [../]

[]

[Kernels]
  [./diff]
    type = MatDiffusion
    variable = C_Mn
    diffusivity = 'diffusion_coefficient'
  [../]
  [./time]
    type = TimeDerivative
    variable = C_Mn
  [../]
[]

[AuxKernels]
  [./ls_oxide_gas]
    type = LineSegmentLevelSetAux
    line_segment_cut_set_user_object = 'moving_cut_oxide'
    variable = ls_oxide_gas
  [../]
  [./ls_metal_oxide]
    type = LineSegmentLevelSetAux
    line_segment_cut_set_user_object = 'fixed_cut_oxide_metal'
    variable = ls_metal_oxide
  [../]
  [./ls_grad]
    type = LineSegmentLevelSetAux
    line_segment_cut_set_user_object = 'fixed_cut_grad'
    variable = ls_grad
  [../]

[]

[Materials]
  [./diffusivity_steel]
    type = GenericConstantMaterial
    prop_names = steel_diffusion_coefficient
    prop_values = 0.36
  [../]
  [./diffusivity_oxide]
    type = GenericConstantMaterial
    prop_names = oxide_diffusion_coefficient
    prop_values = 0.01
  [../]
  [./diffusivity_gas]
    type = GenericConstantMaterial
    prop_names = gas_diffusion_coefficient
    prop_values = 0
  [../]
  [./diff_combined]
    type = LevelSetTriMaterialReal
    levelset_neg_neg_base = 'steel'
    levelset_pos_neg_base = 'oxide'
    levelset_pos_pos_base = 'gas'
    ls_var_1 = ls_metal_oxide
    ls_var_2 = ls_oxide_gas
    prop_name = diffusion_coefficient
    outputs = exodus
  [../]
[]

[BCs]
# Define boundary conditions
  [./left_Mn]
    type = NeumannBC
    variable = C_Mn
    value = 0
    boundary = left
  [../]
  [./left_Mn_init]
    type = DirichletBC
    variable = C_Mn
    value = 7.1445
    boundary = left
  [../]
  [./right_u]
    type = DirichletBC
    variable = C_Mn
    value = 13.3152
    boundary = right
  [../]
[]

[Postprocessors]
  [./position_ox_g]
    type = PositionOfXFEMInterfacePostprocessor
    value_at_interface_uo = value_uo_oxide
    execute_on ='timestep_end final'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type' #-sub_pc_factor_shift_type -pc_factor_shift_amount
  petsc_options_value = 'lu' #NONZERO 1e-6
  line_search = 'none'



  l_tol = 1e-3
  l_max_its = 10
  nl_max_its = 15
  nl_rel_tol = 1e-7
  nl_abs_tol = 1e-7

  start_time = 40
  dt = 10
  num_steps = 46
  max_xfem_update = 1
[]

[Controls]
  [./steady]
    type = TimePeriod
    disable_objects = 'Kernels::time BCs::left_Mn' # Constraints::oxide_metal_constraint
    enable_objects = 'BCs::left_Mn_init UserObjects::value_uo_grad UserObjects::fixed_cut_grad AuxKernels::ls_grad Constraints::grad_constraint' #Constraints::oxide_metal_constraint_init
    start_time = '40'
    end_time = '50'
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
