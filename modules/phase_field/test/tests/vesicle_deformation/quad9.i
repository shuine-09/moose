[GlobalParams]
  epsilon = 0.02
  rz = true
[]

[Problem]
  type = FEProblem
  coord_type = RZ
  rz_coord_axis = y
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 100
  ny = 100
  xmax = 1.0
  ymax = 1.0
  elem_type = QUAD9
[]

[Variables]
  [./u]
#   order = THIRD
#    family = HERMITE
    order = SECOND
    family = LAGRANGE
  [../]
[]

[ICs]
#  [./InitialCondition]
##    type = VesicleIC
##    center = '0.0 0.5 0'
##    variable = u
#    type = SmoothSuperellipsoidIC
#    variable = u
#    x1 = 0.0
#    y1 = 0.375
#    a = 1.0 #0.85  0.35
#    b = 0.25 #0.2   0.85
#    n = 2.0
#    invalue = -1.0
#    outvalue = 1.0
#    int_width = 0.2
#  [../]
[./ic]
  type = VesicleIC
  variable = u
  center = '0.0 0.5 0'
  major = 0.75
  minor = 0.25 #0.3
  epsilon = 0.02
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
    spontaneous_curvature = 0
  [../]
  [./VesiclePenalty]
    type = VesicleVolumeAreaPenalty
    variable = u
    alpha_a = 10000.0 # 500
    alpha_v = 10000.0 # 1000
    vesicle_area = area
    vesicle_volume = volume
    #use_prescribed_area = true
    #prescribed_area = 5.0932644946335
    #use_prescribed_volume = true
    #prescribed_volume = 0.83781841123523
  [../]
[]

[DGKernels]
  [./dg_ch]
    type = DGVesicleShapeDeformation
    variable = u
#   h_elem = 0.01
    eta = 1 #10
  [../]
[]

[BCs]
  [./dirch]
    type = DirichletBC
    variable = u
    boundary = 'bottom right top'
    value = 1.0
  [../]
#  [./flux]
#    type = NeumannBC
#    variable = u
#    boundary = 'left'
#    value = 0.0
#  [../]
#  [./zero_flux]
#    type = DGCHZeroFlux
#    variable = u
#    boundary = 'left'
#    alpha = 10000
#  [../]

[]

[Postprocessors]
  [./volume]
    type = VesicleVolumePostprocessor
    execute_on = 'initial timestep_end'
    variable = u
  [../]
  [./area]
    type = VesicleAreaPostprocessor
    execute_on = 'initial timestep_end'
    variable = u
  [../]
  [./energy]
    type = VesicleTotalEnergyPostprocessor
    execute_on = 'initial timestep_end'
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

  solve_type = 'NEWTON'
#  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
#  petsc_options_value = 'hypre boomeramg 31'

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu       superlu_dist'

#  line_search = 'none'

# petsc_options = '-snes_check_jacobian -sens_check_jacobian_view'

  l_max_its = 8
  l_tol = 1.0e-4
  nl_max_its = 20
  nl_rel_tol = 1.0e-7
  nl_abs_tol = 1.0e-6

  start_time = 0.0
  num_steps = 10000
  dt = 0.00001

#  [./Quadrature]
#     order = second
#  [../]


#  [./TimeStepper]
#    type = SolutionTimeAdaptiveDT
#    dt = 0.0001
#  [../]

#  [./Adaptivity]
#   initial_adaptivity = 1
#    refine_fraction = 0.99
#    coarsen_fraction = 0.01
#    max_h_level = 1
#  [../]

[]

[Outputs]
  print_perf_log = true
  csv = true
  [./exodus]
    type = Exodus
    interval = 20
  [../]
[]
