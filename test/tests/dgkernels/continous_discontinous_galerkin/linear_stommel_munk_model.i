[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 90
  ny = 30
  xmin = 0
  xmax = 3
  ymin = 0
  ymax = 1
  elem_type = QUAD9
  uniform_refine = 1
[]

[Variables]
  active = 'u'
  [./u]
    order =  SECOND
    family = LAGRANGE
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
  active = 'dg_lsm'
  [./dg_lsm]
    type = DGLinearStommelMunk
    variable = u
  [../]
[]

[Functions]
  # A ParsedFunction allows us to supply analytic expressions
  # directly in the input file
  [./exact_sol]
    type = ParsedGradFunction
    value = sin(pi*x/3)^2*sin(pi*y)^2
    grad_x = 2*sin(pi*x/3)*cos(pi*x/3)*pi/3*sin(pi*y)^2
    grad_y = sin(pi*x/3)^2*2*sin(pi*y)*pi*cos(pi*y)
  [../]
[]

[BCs]
  active = 'zero zero_flux'

  [./zero_penalty]
    type = PenaltyDirichletBC
    variable = u
    boundary = '0 1 2 3'
    penalty = 1000
    value = 0
  [../]

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
    alpha = 10
  [../]

[]

[Postprocessors]
  [./dofs]
    type = NumDOFs
  [../]
  [./L2error]
    type = ElementL2Error
    variable = u
    function = exact_sol
  [../]
  [./H1error]
    type = ElementH1Error
    variable = u
    function = exact_sol
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
  solve_type = 'PJFNK'
#  petsc_options = '-snes_mf'
#  petsc_options_iname = '-pc_type -pc_hypre_type'
#  petsc_options_value = 'hypre    boomeramg'

 petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu       superlu_dist'

#  petsc_options = '-snes_check_jacobian -sens_check_jacobian_view'

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
