[GlobalParams]
  order = FIRST
  family = LAGRANGE
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 11
  ny = 11
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1
  elem_type = QUAD4
[]

#[MeshModifiers]
#  [./middle_node]
#    type = AddExtraNodeset
#    new_boundary = 'middle_node'
#    coord = '0.0 0.0'
#  [../]
#[]

[XFEM]
  cut_data = '0.5 1.0 0.5 0.5 0 0'
  qrule = volfrac
  output_cut_plane = true
[]

[Constraints]
  [./xfem_constraint]
    type = XFEMSingleVariableConstraintLS
    variable = u
    level_set_var = ls
    jump = 0
    jump_flux = 0
    diffusivity_ls_plus = 2.0
    diffusivity_ls_minus = 0.5
  [../]
[]

[Variables]
  [./u]
  [../]
[]

[AuxVariables]
  [./ls]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[UserObjects]
  [./xfem_marker_uo]
    type = XFEMMeshCutByLevelSet
    execute_on = 'timestep_end'
    level_set_var = ls
  [../]
[]

[AuxKernels]
  [./ls_function]
    type = XFEMLevelSet
    variable = ls
    function = ls_func
  [../]
[]

[Functions]
  [./u_left]
    type = PiecewiseLinear
    x = '0   2'
    y = '3   5'
  [../]
  [./ls_func]
    type = ParsedFunction
    value = 'x-0.1*(t+1)'
  [../]
[]

[Kernels]
  [./diff]
    type = TwoSidedDiffusion
    variable = u
    diffusivity_ls_plus = 2.0
    diffusivity_ls_minus = 0.01
    level_set = ls
  [../]
[]

[BCs]
# Define boundary conditions
  [./left_u]
    type = DirichletBC
    variable = u
    boundary = 3
    value = 3
  [../]

  [./right_u]
    type = DirichletBC
    variable = u
    boundary = 1
    value = 0
  [../]

#  [./middle_u]
#    type = DirichletBC
#    variable = u
#    boundary = middle_node
#    value = 5
#  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
  line_search = 'none'

  l_tol = 1e-3
  nl_max_its = 15
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-10

  start_time = 0.0
  dt = 1.5
  end_time = 6.0
  max_xfem_update = 1
[]

[Outputs]
  file_base = diffusion_level_set_different_diffusivity_out
  interval = 1
  execute_on = timestep_end
  exodus = true
  [./console]
    type = Console
    perf_log = true
    output_linear = true
  [../]
[]
