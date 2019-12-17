# test for an oxide growing on top of a ziroconium nuclear fuel cladding
# using the C4 model to compute the growth rate

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
  [./velocity]
    type = XFEMC4OxideVelocity
    diffusivity_at_positive_level_set = 0.2978e-5 #1e-5*0.2978
    diffusivity_at_negative_level_set = 0.6667e-11 #1e-11*0.6667
    equilibrium_concentration_jump = 0.3689
    value_at_interface_uo = value_uo
    x0 = 5.9e-4
  [../]
  [./value_uo]
    type = PointValueAtXFEMInterface
    variable = 'u'
    geometric_cut_userobject = 'moving_line_segments'
    execute_on = 'nonlinear'
    level_set_var = ls
  [../]
  [./moving_line_segments]
    type = MovingLineSegmentCutSetUserObject
    cut_data = '5.9e-4 0 5.9e-4 1e-3 0 0'
    heal_always = true
    interface_velocity = velocity
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
    function = 'if(x<5.9e-4, 0.03, 0.6667)'
  [../]
[]

[AuxVariables]
  [./ls]
    order = FIRST
    family = LAGRANGE
  [../]
[]

# Need to use XFEMTwoSideDirichlet that has been removed
[Constraints]
  [./u_constraint]
    type = XFEMTwoSideDirichlet
    geometric_cut_userobject = 'moving_line_segments'
    use_displaced_mesh = false
    variable = u
    value_at_positive_level_set_interface = 1 #0.2978
    value_at_negative_level_set_interface = 1 #0.6667
#    value = 0.6667
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
  [./ls]
    type = LineSegmentLevelSetAux
    line_segment_cut_set_user_object = 'moving_line_segments'
    variable = ls
  [../]
[]

[Materials]
  [./diffusivity_A]
    type = GenericConstantMaterial
    prop_names = A_diffusion_coefficient
    prop_values = 1e-5
  [../]
  [./diffusivity_B]
    type = GenericConstantMaterial
    prop_names = B_diffusion_coefficient
    prop_values = 1e-11
  [../]
  [./diff_combined]
    type = LevelSetBiMaterialReal
    levelset_positive_base = 'A'
    levelset_negative_base = 'B'
    level_set_var = ls
    prop_name = diffusion_coefficient
  [../]
[]

[BCs]
# Define boundary conditions
  [./left_u]
    type = DirichletBC
    variable = u
    value = 0.0450 #0.03/0.6667
    boundary = left
  [../]

  [./right_u]
    type = DirichletBC
    variable = u
    value = 2.2388 #0.6667/0.2978
    boundary = right
  [../]
[]

[Postprocessors]
  [./cut_data_x]
    type = OxideLayerThickness
    moving_line_segments = moving_line_segments
    execute_on = TIMESTEP_END
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
