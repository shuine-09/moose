[GlobalParams]
  order = FIRST
  family = LAGRANGE
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 21
  ny = 11
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
    cut_data = '0.5 0 0.5 0.5 0 0
                0.5 0.5 0.5 1.0 0 0'
    # cut_data = '0.3 0 0.5 0.5 0 0
    #             0.5 0.5 0.7 1.0 0 0'
    # cut_data = '  1.4999    1.0087    1.4908    1.0954         0         0
    # 1.4908    1.0954    1.4668    1.1792         0         0
    # 1.4668    1.1792    1.4286    1.2575         0         0
    # 1.4286    1.2575    1.3774    1.3280         0         0
    # 1.3774    1.3280    1.3147    1.3886         0         0
    # 1.3147    1.3886    1.2424    1.4373         0         0
    # 1.2424    1.4373    1.1628    1.4728         0         0
    # 1.1628    1.4728    1.0782    1.4938         0         0
    # 1.0782    1.4938    0.9913    1.4999         0         0
    # 0.9913    1.4999    0.9046    1.4908         0         0
    # 0.9046    1.4908    0.8208    1.4668         0         0
    # 0.8208    1.4668    0.7425    1.4286         0         0
    # 0.7425    1.4286    0.6720    1.3774         0         0
    # 0.6720    1.3774    0.6114    1.3147         0         0
    # 0.6114    1.3147    0.5627    1.2424         0         0
    # 0.5627    1.2424    0.5272    1.1628         0         0
    # 0.5272    1.1628    0.5062    1.0782         0         0
    # 0.5062    1.0782    0.5001    0.9913         0         0
    # 0.5001    0.9913    0.5092    0.9046         0         0
    # 0.5092    0.9046    0.5332    0.8208         0         0
    # 0.5332    0.8208    0.5714    0.7425         0         0
    # 0.5714    0.7425    0.6226    0.6720         0         0
    # 0.6226    0.6720    0.6853    0.6114         0         0
    # 0.6853    0.6114    0.7576    0.5627         0         0
    # 0.7576    0.5627    0.8372    0.5272         0         0
    # 0.8372    0.5272    0.9218    0.5062         0         0
    # 0.9218    0.5062    1.0087    0.5001         0         0
    # 1.0087    0.5001    1.0954    0.5092         0         0
    # 1.0954    0.5092    1.1792    0.5332         0         0
    # 1.1792    0.5332    1.2575    0.5714         0         0
    # 1.2575    0.5714    1.3280    0.6226         0         0
    # 1.3280    0.6226    1.3886    0.6853         0         0
    # 1.3886    0.6853    1.4373    0.7576         0         0
    # 1.4373    0.7576    1.4728    0.8372         0         0
    # 1.4728    0.8372    1.4938    0.9218         0         0
    # 1.4938    0.9218    1.4999    1.0087         0         0'
    heal_always = true
    var = u
    interface_value_uo = value_uo
  [../]
[]

[Variables]
  [./u]
    initial_condition = 10
  [../]
[]

[AuxVariables]
  [./ls]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Constraints]
  [./u_constraint]
    type = XFEMTwoSideDirichlet
    geometric_cut_userobject = 'moving_line_segments'
    use_displaced_mesh = false
    variable = u
    level_set_positive_value = 10
    level_set_negative_value = 2
    #level_set_negative_value = 10
    alpha = 1e6
  [../]
[]

[Kernels]
  [./diff]
    type = Diffusion
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
[]

[BCs]
# Define boundary conditions
  [./left_u]
    type = DirichletBC
    variable = u
    value = 10
    boundary = 3
  [../]

  # [./left_u_all]
  #   type = XFEMOneSideConstantBC
  #   variable = u
  #   value = 10
  #   left_x = 0.47
  #   boundary = all
  # [../]

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
  dt = 0.05
  end_time = 1.5
  num_steps = 3
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
