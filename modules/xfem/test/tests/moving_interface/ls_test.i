[GlobalParams]
  order = FIRST
  family = LAGRANGE
[]

[XFEM]
  geometric_cut_userobjects = 'cut_mesh'
  qrule = volfrac
  output_cut_plane = true
[]

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 50#15
    ny = 50#15
    nz = 50#15
    xmin = -2
    xmax = 2
    ymin = -2
    ymax = 2
    zmin = -2
    zmax = 2
    elem_type = HEX8
  []
[]

[UserObjects]
  [./cut_mesh]
    type = InterfaceMeshCut3DUserObject
    mesh_file = bunny.xda
    velocity = 0.002
    heal_always = true
  [../]
[]

# [UserObjects]
#   [./level_set_cut_uo]
#     type = LevelSetCutUserObject
#     level_set_var = ls
#   [../]
# []
#
# [AuxVariables]
#   [./ls]
#     order = FIRST
#     family = LAGRANGE
#   [../]
# []
#
# [AuxKernels]
#   [./ls_function]
#     type = FunctionAux
#     variable = ls
#     function = ls_func
#   [../]
# []

[Functions]
  [./ls_func]
    type = ParsedFunction
    value = 'sqrt(x*x+y*y+z*z)-0.4'
  [../]
[]

[Variables]
  [./u]
  [../]
[]

[AuxVariables]
  [./ls]
  [../]
[]

[AuxKernels]
  [./ls]
    type = MeshCutLevelSetAux
    mesh_cut_user_object = cut_mesh
    variable = ls
    execute_on = 'TIMESTEP_END'
  [../]
[]

[Kernels]
  [./diff]
    type = MatDiffusion
    variable = u
    diffusivity = 1
  [../]
  [./time_deriv]
    type = TimeDerivative
    variable = u
  [../]
[]

[BCs]
  [./front_u]
    type = DirichletBC
    variable = u
    boundary = front
    value = 0
  [../]
  [./back_u]
    type = DirichletBC
    variable = u
    boundary = back
    value = 1
  [../]
[]

[Executioner]
  type = Transient

  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  l_max_its = 20
  l_tol = 1e-3
  nl_max_its = 15
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-10

  start_time = 0.0
  dt = 1
  end_time = 10

  max_xfem_update = 1
[]

[Outputs]
  exodus = true
  execute_on = timestep_end
  csv = true
  perf_graph = true
  [./console]
    type = Console
    output_linear = true
  [../]
[]
