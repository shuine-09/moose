[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 121
  ny = 25
  xmin = 0
  xmax = 30
  ymin = 0
  ymax = 6
  elem_type = QUAD9
[]

[XFEM]
  geometric_cut_userobjects = 'level_set_uo'
  qrule = volfrac
  output_cut_plane = true
[]

[AuxVariables]
  [./ls]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxKernels]
  [./ls_function]
    type = FunctionAux
    variable = ls
    #function = 'sqrt((x-15)*(x-15)+(y-3)*(y-3))-1.9'
    #function = 'sqrt((x-10)*(x-10)+(y-3)*(y-3))-1.9'
    function = 'sqrt((x-10)*(x-10)/1.31/1.31+(y-3)*(y-3)/1.91/1.91)-1'
  [../]
[]

[UserObjects]
  [./level_set_uo]
    type = LevelSetCutUserObject
    level_set_var = ls
    time_start_cut = 0.0
    time_end_cut = 0.0
  [../]
[]

[GlobalParams]
  gravity = '0 -0.001 0'
  convective_term = false
  integrate_p_by_parts = false
[]

[MeshModifiers]
  [./corner_node]
    type = AddExtraNodeset
    new_boundary = top_right
    coord = '30 6'
  [../]
[]

[Variables]
  [./vel_x]
    order = SECOND
    family = LAGRANGE
  [../]
  [./vel_y]
    order = SECOND
    family = LAGRANGE
  [../]
  [./p]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./mass]
    type = INSMass
    variable = p
    u = vel_x
    v = vel_y
  [../]
  [./x_momentum_space]
    type = INSMomentumLaplaceForm
    variable = vel_x
    u = vel_x
    v = vel_y
    p = p
    component = 0
  [../]
  [./y_momentum_space]
    type = INSMomentumLaplaceForm
    variable = vel_y
    u = vel_x
    v = vel_y
    p = p
    component = 1
  [../]
[]

[BCs]
  [./x_no_slip]
    type = DirichletBC
    variable = vel_x
    boundary = 'top bottom'
    value = 0.0
  [../]
  [./y_no_slip]
    type = DirichletBC
    variable = vel_y
    boundary = 'left top bottom'
    value = 0.0
  [../]
  [./x_inlet]
    type = FunctionDirichletBC
    variable = vel_x
    boundary = 'left'
    function = 'inlet_func'
  [../]
  [./p_corner]
    # Since the pressure is not integrated by parts in this example,
    # it is only specified up to a constant by the natural outflow BC.
    # Therefore, we need to pin its value at a single location.
    type = DirichletBC
    boundary = top_right
    value = 0
    variable = p
  [../]
[]

[Constraints]
  [./xfem_constraint_x]
    type = XFEMTwoSideDirichletBC
    variable = vel_x
    alpha = 1000
  [../]
  [./xfem_constraint_y]
    type = XFEMTwoSideDirichletBC
    variable = vel_y
    alpha = 1000
  [../]
[]

[Materials]
  [./const]
    type = GenericConstantMaterial
    prop_names = 'rho mu'
    prop_values = '100  1'
  [../]
[]

[Preconditioning]
  [./SMP_PJFNK]
    type = SMP
    full = true
    solve_type = NEWTON
  [../]
[]

[VectorPostprocessors]
  [./x7a]
    type = LineValueSampler
    variable = 'p vel_x vel_y'
    start_point = '7 0 0'
    end_point = '7 0.999 0'
    num_points = '20'
    sort_by = y
  [../]
  [./x7b]
    type = LineValueSampler
    variable = 'p vel_x vel_y'
    start_point = '7 4.001 0'
    end_point = '7 6 0'
    num_points = '40'
    sort_by = y
  [../]
  [./y5]
    type = LineValueSampler
    variable = 'p vel_x vel_y'
    start_point = '0 5 0'
    end_point = '30 5 0'
    num_points = '500'
    sort_by = x
  [../]
[]

[Executioner]
  type = Transient
#  petsc_options_iname = '-ksp_gmres_restart -pc_type -sub_pc_type -sub_pc_factor_levels'
#  petsc_options_value = '300                bjacobi  ilu          4'
petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
petsc_options_value = 'lu mumps'
  line_search = none
  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-10
  nl_max_its = 15
  l_tol = 1e-10
  l_max_its = 1000
  dt = 1.0;
  end_time = 1.0;
[]

[Outputs]
  exodus = true
  csv = true
  execute_on = TIMESTEP_END
[]

[Functions]
  [./inlet_func]
    type = ParsedFunction
    value = '-0.001 * (y-3)^2 + 0.009'
  [../]
[]
