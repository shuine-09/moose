[GlobalParams]
  epsilon = 0.04
[]

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 30
  ny = 30
  nz = 60
  xmax = 0.6
  ymax = 0.6
  zmax = 1.2
  elem_type = HEX27
[]

[Variables]
  [./u]
    order = SECOND
    family = LAGRANGE
  [../]
[]

[ICs]
  [./InitialCondition]
    type = VesicleIC
    center = '0.3 0.3 0.6'
    variable = u
#    type = SmoothSuperellipsoidIC
#    variable = u
#    x1 = 0.5
#    y1 = 0.5
#    z1 = 0.5
#    a = 0.2
#    b = 0.4
#    n = 2.0
#    invalue = -1.0
#    outvalue = 1.0
#    int_width = 0.1
  [../]
[]

[Kernels]
  [./ie_c]
    type = TimeDerivative
    variable = u
  [../]
  [./VesicleShape]
    type = VesicleShapeDeformation
    variable = u
  [../]
  [./VesiclePenalty]
    type = VesicleVolumeAreaPenalty
    variable = u
    alpha_a = 1000.0
    alpha_v = 10000.0
    vesicle_area = area
    vesicle_volume = volume
  [../]
[]

[DGKernels]
  [./dg_ch]
    type = DGVesicleShapeDeformation
    variable = u
    eta = 100
  [../]
[]

[BCs]
  [./dirch]
    type = DirichletBC
    variable = u
    boundary = 'back front left right bottom top'
    value = 1.0
  [../]
#  [./flux]
#    type = NeumannBC
#    variable = u
#    boundary = 'bottom right top'
#    value = 0.0
#  [../]
#  [./zero_flux]
#    type = DGCHZeroFlux
#    variable = u
#    boundary = 'bottom right top'
#    alpha = 100
#  [../]

[]

[Postprocessors]
  [./volume]
    type = VesicleVolumePostprocessor
    execute_on = 'initial nonlinear timestep_end'
    variable = u
  [../]
  [./area]
    type = VesicleAreaPostprocessor
    execute_on = 'initial nonlinear timestep_end'
    variable = u
  [../]
[]


[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
#  scheme = 'bdf2'

  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 31'

#  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
#  petsc_options_value = 'lu       superlu_dist'


#  petsc_options = '-snes_check_jacobian -sens_check_jacobian_view'

  l_max_its = 15
  l_tol = 1.0e-4
  nl_max_its = 20
  nl_rel_tol = 1.0e-7
  nl_abs_tol = 1.0e-5

  start_time = 0.0
  num_steps = 2000
  dt = 0.00001

#  [./Adaptivity]
#   initial_adaptivity = 1
#    refine_fraction = 0.4
#    coarsen_fraction = 0.2
#    max_h_level = 1
#  [../]

[]

[Outputs]
  exodus = true
  csv = true
[]
