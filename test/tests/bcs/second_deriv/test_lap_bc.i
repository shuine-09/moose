[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 1
  ny = 1
  xmin = -1
  xmax = 1
  ymin = -1
  ymax = 1
  elem_type = QUAD9
[]

[Variables]
  active = 'u'

  [./u]
    order = SECOND
    family = LAGRANGE
  [../]
[]

[Kernels]
  active = 'diff'

  [./diff]
    type = Diffusion
    variable = u
  [../]
  [./ffn]
    type = UserForcingFunction
    variable = u
    function = force_fn
  [../]
[]

[Functions]
  [./left_bc_func]
    type = ParsedFunction
    value = '1+y*y'
  [../]
  [./right_bc_func]
    type = ParsedFunction
    value = '1+y*y'
  [../]
  [./top_bc_func]
    type = ParsedFunction
    value = '1+x*x'
  [../]
  [./bottom_bc_func]
    type = ParsedFunction
    value = '1+x*x'
  [../]
  [./force_fn]
    type = ParsedFunction
    value = -4
  [../]
[] 

[BCs]
  active = 'left top bottom right_test'
  [./left]
    type = FunctionDirichletBC
    variable = u
    boundary = left
    function = left_bc_func
  [../]
  [./right]
    type = FunctionDirichletBC
    variable = u
    boundary = right
    function = right_bc_func
  [../]
  [./bottom]
    type = FunctionDirichletBC
    variable = u
    boundary = bottom
    function = bottom_bc_func
  [../]
  [./top]
    type = FunctionDirichletBC
    variable = u
    boundary = top
    function = top_bc_func
  [../]
  [./right_test]
    type = TestLapBC
    variable = u
    boundary = right
  [../]
[]

[Executioner]
  type = Steady

  solve_type = 'NEWTON'

  # Force 2x2 quadrature on QUAD9 just to cut down the number of values
  # there are to look at.
  [./Quadrature]
    type = GAUSS
    order = SECOND
  [../]

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  petsc_options = '-snes_check_jacobian -snes_check_jacobian_view'
[]

[Outputs]
  file_base = out
  exodus = true
[]
