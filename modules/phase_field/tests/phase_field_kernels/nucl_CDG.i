#
# Test the non-split parsed function free enery Cahn-Hilliard Bulk kernel
# The free energy used here has the same functional form as the CHPoly kernel
# If everything works, the output of this test should replicate the output
# of marmot/tests/chpoly_test/CHPoly_test.i (exodiff match)
#

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 128
  ny = 128
  xmin = 0
  xmax = 256
  ymin = 0
  ymax = 256
  elem_type = QUAD9
#  uniform_refine = 2
[]

[Variables]
  [./cv]
    order = SECOND
    family = LAGRANGE
  [../]
[]

[ICs]
#  [./InitialCondition]
#    type = CrossIC
#    x1 = 5.0
#    y1 = 5.0
#    x2 = 45.0
#    y2 = 45.0
#    type = RandomIC
#    min = -0.6
#    max = -0.4
#    variable = cv
  [./conIC_c]
    type = ConstantIC
    variable = cv
    value = 0.1
  [../]
[]

[Kernels]
  [./ie_c]
    type = TimeDerivative
    variable = cv
  [../]
  [./CHSolid]
    type = CahnHilliard
    variable = cv
    f_name = F
    mob_name = M
  [../]
  [./CHInterface]
    type = CHInterface
    variable = cv
    mob_name = M
    kappa_name = kappa_c
  [../]
[]

[DGKernels]
  active =  'dg_ch'
  [./dg_ch]
    type = DGCahnHilliard
    variable = cv
    eta = 50
    kappa_name = kappa_c
    mob_name = M
  [../]
[]

[Materials]
#  [./consts]
#    type = GenericConstantMaterial
#    prop_names  = 'M kappa_c'
#    prop_values = '1 0.5'
#    block = 0
#  [../]
#  [./free_energy]
#    type = DerivativeParsedMaterial
#    block = 0
#    f_name = F
#    args = 'cv'
##    function = '(1-cv)^2 * (1+cv)^2'
#    function = 'cv^4/12 - cv^2/2'
#  [../]
  [./FreeEng]
      type = DerivativeParsedMaterial
      function = '(A*(cv-c1)^2*(cv-c2)^2)'
      args = 'cv'
      f_name = F
      constant_names = 'A c1 c2'
      constant_expressions = '39.0625 0.1 0.9'
      derivative_order = 2
    [../]
  [./const]
    type = GenericConstantMaterial
    prop_names = 'M kappa_c'
    prop_values = '0.1 64'
  [../]
[]

[BCs]
#  [./Periodic]
#    [./all]
#       auto_direction = 'x y'
#    [../]
#  [../]
  
  [./zero_flux]
    type = DGCHZeroFlux
    variable = cv
    boundary = '0 1 2 3'
    alpha = 100
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[DiracKernels]
  [./radiation]
    type = RandomPointSource
    num_srcs = 256
     value = 0.25
    variable = cv
    x1 = 0.0
    y1 = 0.0
    x2 = 256.0
    y2 = 256.0
  [../]
[]


[Executioner]
  type = Transient
  scheme = 'bdf2'

  solve_type = 'NEWTON'
#  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
#  petsc_options_value = 'hypre boomeramg 31'

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu       superlu_dist'

#  petsc_options = '-snes_check_jacobian -sens_check_jacobian_view'

  l_max_its = 15
  l_tol = 1.0e-4
  nl_max_its = 20
  nl_rel_tol = 1.0e-8
  nl_abs_tol = 1.0e-8

  start_time = 0.0
  num_steps = 1000
  [./TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = .75
    dt = 0.5
    growth_factor = 1.2
    optimal_iterations = 8
  [../]
#  [./Adaptivity]
#    refine_fraction = 0.5
#    coarsen_fraction = 0.0
#    max_h_level = 2
#    initial_adaptivity = 2
#  [../]
[]

[Outputs]
  [./out]
    type = Exodus
# refinements = 1
  [../]
[]
