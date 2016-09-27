[GlobalParams]
  order = FIRST
  family = LAGRANGE
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 81
  ny = 81
  xmin = 0
  xmax = 1.0
  ymin = 0
  ymax = 1.0
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
  interface_points = '0.05 1.0 0
                      0.05 0.9 0
                      0.05 0.8 0
                      0.05 0.7 0
                      0.05 0.6 0
                      0.05 0.5 0
                      0.05 0.4 0
                      0.05 0.3 0
                      0.05 0.2 0
                      0.05 0.1 0
                      0.05 0   0'
  qrule = volfrac
  output_cut_plane = true
[]

[Constraints]
  [./xfem_constraint]
    type = XFEMSingleVariableConstraintPts
    variable = u
    jump = 0
    jump_flux = 0
  [../]
[]

[UserObjects]
  [./xfem_moving_points]
    type = XFEMMovingInterfacePoints
    execute_on = timestep_end
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

[Functions]
  [./u_left]
    type = PiecewiseLinear
    x = '0   2'
    y = '3   5'
  [../]
  [./u_left2]
    type = ParsedFunction
    value = '-2 * (y-0.5) * (y-0.5) + 1.0'
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
#  [./left_u]
#    type = DirichletBC
#    variable = u
#    boundary = 3
#    value = 3
#  [../]
  
  [./left_u2]
    type = FunctionDirichletBC
    variable = u
    boundary = 3
    function = u_left2
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
#  petsc_options_iname = '-pc_type -pc_hypre_type'
#  petsc_options_value = 'hypre boomeramg'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  line_search = 'none'

  l_tol = 1e-3
  nl_max_its = 15
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-10

  start_time = 0.0
  dt = 1.0
  end_time = 30
  max_xfem_update = 1
[]

[Outputs]
  file_base = diffusion_cut_points_out2
  interval = 1
  execute_on = timestep_end
  exodus = true
  [./console]
    type = Console
    perf_log = true
    output_linear = true
  [../]
[]
