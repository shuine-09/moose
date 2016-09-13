[GlobalParams]
  order = FIRST
  family = LAGRANGE
[]

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 41
  ny = 41
  nz = 41
  xmin = 0.0
  xmax = 1.0
  ymin = 0.0
  ymax = 1.0
  zmin = 0.0
  zmax = 1.0
  elem_type = HEX8
[]

[XFEM]
  cut_data = '0.5 1.0 0.5 0.5 0 0'
  qrule = volfrac
  output_cut_plane = true
[]

[Variables]
  [./u]
  [../]
[]

[Functions]
  [./u_left]
    type = PiecewiseLinear
    x = '0   2'
    y = '0  0.1'
  [../]
  [./ls_func]
    type = ParsedFunction
    value = 'sqrt(x*x + y*y + z*z) - 0.1*(t+1.1)'
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

[Kernels]
  [./diff]
    type = Diffusion
    variable = u
  [../]
[]

[BCs]
# Define boundary conditions
  [./left_u]
    type = DirichletBC
    variable = u
    boundary = left
    value = 2
  [../]

  [./right_u]
    type = DirichletBC
    variable = u
    boundary = right
    value = 1
  [../]
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
  dt = 0.25
  end_time = 10.0
  max_xfem_update = 1

[]

[Outputs]
  file_base = diffusion_level_set_3d_out
  interval = 1
  execute_on = timestep_end
  exodus = true
  [./console]
    type = Console
    perf_log = true
    output_linear = true
  [../]
[]
