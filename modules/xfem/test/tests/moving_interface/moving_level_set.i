[GlobalParams]
  order = FIRST
  family = LAGRANGE
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
<<<<<<< HEAD:modules/xfem/test/tests/moving_interface/moving_level_set.i
  nx = 5
  ny = 5
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1
  elem_type = QUAD4
[]

[XFEM]
  qrule = volfrac
  output_cut_plane = true
[]

[UserObjects]
  [./line_seg_cut_uo]
    type = LineSegmentCutSetUserObject
    cut_data = '0.3 1.0 0.3 0.2 0 3'
    heal_always = false
  [../]
  [./level_set_cut_uo]
    type = LevelSetCutUserObject
    level_set_var = ls
    heal_always = true
  [../]
[]

[Variables]
  [./u]
  [../]
=======
  nx = 11
  ny = 1
  xmin = 0
  xmax = 2
  ymin = 0
  ymax = 1
  elem_type = QUAD4
>>>>>>> preliminary development of moving interface capability using XFEM:modules/xfem/test/tests/moving_interface/phase_transition.i
[]

[XFEM]
  qrule = volfrac
  output_cut_plane = true
[]

[UserObjects]
  [./velocity]
    type = XFEMPhaseTransitionMovingInterfaceVelocity
    diffusivity_at_positive_level_set = 5
    diffusivity_at_negative_level_set = 1
    equilibrium_concentration_jump = 1
    value_at_interface_uo = value_uo
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
    cut_data = '0.5 0 0.5 1.0 0 0'
    heal_always = true
    interface_velocity = velocity
  [../]
[]

<<<<<<< HEAD:modules/xfem/test/tests/moving_interface/moving_level_set.i
[Functions]
  [./u_left]
    type = PiecewiseLinear
    x = '0   2'
    y = '3   5'
  [../]
  [./ls_func]
    type = ParsedFunction
    value = 'x-0.7-0.07*(t-1)'
=======
[Variables]
  [./u]
  [../]
[]

[ICs]
  [./ic_u]
    type = FunctionIC
    variable = u
    function = 'if(x<0.51, 2, 1)'
  [../]
[]

[AuxVariables]
  [./ls]
    order = FIRST
    family = LAGRANGE
>>>>>>> preliminary development of moving interface capability using XFEM:modules/xfem/test/tests/moving_interface/phase_transition.i
  [../]
[]

[Constraints]
  [./u_constraint]
<<<<<<< HEAD:modules/xfem/test/tests/moving_interface/moving_level_set.i
    type = XFEMSingleVariableConstraint
    geometric_cut_userobject = 'level_set_cut_uo'
    use_displaced_mesh = false
    variable = u
    use_penalty = true
=======
    type = XFEMEqualValueAtInterface
    geometric_cut_userobject = 'moving_line_segments'
    use_displaced_mesh = false
    variable = u
    value = 2
>>>>>>> preliminary development of moving interface capability using XFEM:modules/xfem/test/tests/moving_interface/phase_transition.i
    alpha = 1e5
  [../]
[]

[Kernels]
  [./diff]
<<<<<<< HEAD:modules/xfem/test/tests/moving_interface/moving_level_set.i
    type = Diffusion
    variable = u
  [../]
[]

[BCs]
# Define boundary conditions
  [./left_u]
    type = DirichletBC
    variable = u
    boundary = 3
    value = 3
=======
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
>>>>>>> preliminary development of moving interface capability using XFEM:modules/xfem/test/tests/moving_interface/phase_transition.i
  [../]

<<<<<<< HEAD:modules/xfem/test/tests/moving_interface/moving_level_set.i
  [./right_u]
    type = DirichletBC
    variable = u
    boundary = 1
    value = 0
  [../]

=======
[Materials]
  [./diffusivity_A]
    type = GenericConstantMaterial
    prop_names = A_diffusion_coefficient
    prop_values = 5
  [../]
  [./diffusivity_B]
    type = GenericConstantMaterial
    prop_names = B_diffusion_coefficient
    prop_values = 1
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
    value = 2
    boundary = 3
  [../]

  [./right_u]
    type = NeumannBC
    variable = u
    boundary = 1
    value = 0
  [../]
>>>>>>> preliminary development of moving interface capability using XFEM:modules/xfem/test/tests/moving_interface/phase_transition.i
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
<<<<<<< HEAD:modules/xfem/test/tests/moving_interface/moving_level_set.i
  # petsc_options_iname = '-pc_type -pc_hypre_type'
  # petsc_options_value = 'hypre boomeramg'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
=======
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
>>>>>>> preliminary development of moving interface capability using XFEM:modules/xfem/test/tests/moving_interface/phase_transition.i
  line_search = 'none'

  l_tol = 1e-3
  nl_max_its = 15
<<<<<<< HEAD:modules/xfem/test/tests/moving_interface/moving_level_set.i
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-9

  start_time = 0.0
  dt = 1
  end_time = 3.0
=======
  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-11

  start_time = 0.0
  dt = 0.01
  num_steps = 4
>>>>>>> preliminary development of moving interface capability using XFEM:modules/xfem/test/tests/moving_interface/phase_transition.i
  max_xfem_update = 1
[]


[Outputs]
<<<<<<< HEAD:modules/xfem/test/tests/moving_interface/moving_level_set.i
  interval = 1
=======
>>>>>>> preliminary development of moving interface capability using XFEM:modules/xfem/test/tests/moving_interface/phase_transition.i
  execute_on = timestep_end
  exodus = true
  [./console]
    type = Console
    output_linear = true
  [../]
  csv = true
[]
