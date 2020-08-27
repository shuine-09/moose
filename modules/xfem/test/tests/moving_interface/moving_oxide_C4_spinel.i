# Input file for an oxide growing outward on top of a steel 21-2N sample
# using the C4 model to compute the growth rate
# The variable is the Mn concentration [/nm^3]
# The length unit is the nanometer. The time unit is the hour
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
  nx = 41
  ny = 20
  xmin = 0
  xmax = 8000
  ymin = 0
  ymax = 4000
  elem_type = QUAD4
[]

[XFEM]
  qrule = volfrac
  output_cut_plane = true
[]

[UserObjects]
  [./fixed_cut_oxide_metal]
    type = LineSegmentCutSetUserObject
    cut_data = '5000 0 5000 4000 0 0'
  [../]
  [./velocity_oxide]
    type = XFEMC4VelocitySteelOx
    value_at_interface_uo = value_uo_oxide
  [../]
  [./value_uo_oxide]
    type = PointValueAtXFEMInterface
    variable = 'C_Mn'
    geometric_cut_userobject = 'moving_cut_oxide'
    execute_on = 'nonlinear'
    level_set_var = ls_oxide_gas
  [../]
  [./moving_cut_oxide]
    type = MovingLineSegmentCutSetUserObject
    cut_data = '5500 0 5500 4000 0 0'
    heal_always = true
    interface_velocity = velocity_oxide
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
    function = 'if(x<5000, 7.1445,if(x<5500,13.3179,13.3152))'
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
[]


[Constraints]
  [./oxide_metal_constraint]
    type = XFEMTwoSideDirichlet
    geometric_cut_userobject = 'fixed_cut_oxide_metal'
    use_displaced_mesh = false
    variable = C_Mn
    value_at_positive_level_set_interface = 2.1577
    value_at_negative_level_set_interface = 13.3179
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
[]

[Materials]
  [./diffusivity_steel]
    type = GenericConstantMaterial
    prop_names = steel_diffusion_coefficient
    prop_values = 36000
  [../]
  [./diffusivity_oxide]
    type = GenericConstantMaterial
    prop_names = oxide_diffusion_coefficient
    prop_values = 10000
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

  start_time = 40
  dt = 10
  num_steps = 46
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
