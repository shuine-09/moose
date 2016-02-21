[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 30
  ny = 10
  xmin = 0
  xmax = 3
  ymin = 0
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

[Functions]
  active = 'forcing_fn'

  [./forcing_fn]
    type = ParsedFunction
#    function = -4.0+(x*x)+(y*y)
#    function = x
#    function = (x*x)-2.0
#    value = 2*pow(e,-x-(y*y))*(1-2*y*y)
     value = '-0.2-2*x*x -10*sin(pi*y)'
#    function = (x*x*x)-6.0*x
  [../]
[]

[Kernels]
  active = 'lsm forcing'
  [./lsm]
    type = LinearStommelMunk
    variable = u
  [../]

  [./forcing]
    type = LSMForcingFunction
    variable = u
  [../]
[]

[DGKernels]
  active = ''
  [./dg_lsm]
    type = DGLinearStommelMunk
    variable = u
  [../]
[]

[BCs]
  active = 'zero zero_flux'

  [./zero]
    type = DirichletBC
    variable = u
    boundary = '0 1 2 3'
    value = 0
  [../]

  [./zero_flux]
    type = LSMDirichletBC
    variable = u
    boundary = '0 1 2 3'
    alpha = 1
  [../]

[]


[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Steady

  # Preconditioned JFNK (default)
  solve_type = 'NEWTON'
#  petsc_options = '-snes_mf'
#  petsc_options_iname = '-pc_type -pc_hypre_type'
#  petsc_options_value = 'hypre    boomeramg'

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu       superlu_dist'

  petsc_options = '-snes_check_jacobian -sens_check_jacobian_view'

#  petsc_options = '-snes_mf'
#  max_r_steps = 2
 
   nl_rel_tol = 1e-5
  
   l_tol = 1e-5

#  nl_rel_tol = 1e-12

   [./Quadrature]
     order = second
   [../]
[]

[Outputs]
  file_base = out
  exodus = true
  csv = true
[]
