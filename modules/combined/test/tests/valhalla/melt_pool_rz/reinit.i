[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 0.005
    ymin = 0
    ymax = 0.01
    nx = 200
    ny = 400
    elem_type = QUAD4
  []
[]

[Problem]
  coord_type = RZ
[]

# [Adaptivity]
#   steps = 1
#   marker = box
#   max_h_level = 1
#   initial_steps = 1
#   stop_time = 1.0e-8
#   [Markers]
#     [box]
#       bottom_left = '0.000 0.004 0'
#       inside = refine
#       top_right = '0.0025 0.006 0'
#       outside = do_nothing
#       type = BoxMarker
#     []
#   []
# []

[Variables]
  [ls]
    order = FIRST
  []
  # [grad_ls]
  #   family = LAGRANGE_VEC
  # []
[]

[AuxVariables]
  [ls_0]
    order = FIRST
  []
[]

[Kernels]
  # [grad_ls]
  #   type = VariableGradientRegularization
  #   regularized_var = ls_0
  #   variable = grad_ls
  # []
  [time]
    type = TimeDerivative
    variable = ls
  []
  # [reinit]
  #   type = LevelSetGradientRegularizationReinitialization
  #   variable = ls
  #   level_set_gradient = grad_ls
  #   epsilon = 0.00004
  # []
  [reinit]
    type = LevelSetOlssonReinitialization
    variable = ls
    phi_0 = ls_0
    epsilon = 0.00004
  []
[]

[Problem]
  type = LevelSetReinitializationProblem
[]

[UserObjects]
  [arnold]
    type = LevelSetOlssonTerminator
    tol = 0.5
    min_steps = 5
  []
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
    solve_type = 'NEWTON'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  start_time = 0
  num_steps = 10
  nl_abs_tol = 1e-14
  nl_max_its = 10
  line_search = none
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu superlu_dist'
  dt = 0.00001
[]

[Outputs]
  exodus = false
  execute_on = 'TIMESTEP_END'
[]
