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
[]

[BCs]
  [./left]
    type = DirichletBC
    variable = u
    boundary = left
    value = 1
  [../]

  [./right]
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
