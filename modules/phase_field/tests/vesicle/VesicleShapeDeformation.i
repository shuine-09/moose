[GlobalParams]
  epsilon = 0.03;
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 200
  ny = 100
  xmax = 2
  ymax = 1
  elem_type = QUAD9
[]

[Variables]
  [./u]
    order = SECOND
    family = LAGRANGE
  [../]
[]

[ICs]
  [./InitialCondition]
#    type = VesicleIC
#    center = '0.5 0.5 0'
#    variable = u
    type = SmoothSuperellipsoidIC
    variable = u
    x1 = 1.0
    y1 = 0.5
    a = 0.7
    b = 0.25
    n = 2.5
    invalue = -1.0
    outvalue = 1.0
    int_width = 0.1
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

#[UserObjects]
#  [./vesicle_volume]
#    type = VesicleVolume
#    variable = u
#    mesh_volume = Volume
#    execute_on = 'initial timestep_end'
#  [../]
#  [./vesicle_area]
#    type = VesicleArea
#    variable = u
#    execute_on = 'initial timestep_end'
#  [../]
#[]

[BCs]
  [./dirch]
    type = DirichletBC
    variable = u
    boundary = '0 1 2 3'
    value = 1
  [../]
#  [./zero_flux]
#    type = DGCHZeroFlux
#    variable = u
#    boundary = '0 1 2 3'
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
#  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
#  petsc_options_value = 'hypre boomeramg 31'

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu       superlu_dist'


#  petsc_options = '-snes_check_jacobian -sens_check_jacobian_view'

  l_max_its = 15
  l_tol = 1.0e-4
  nl_max_its = 20
  nl_rel_tol = 1.0e-7
  nl_abs_tol = 1.0e-5

  start_time = 0.0
  num_steps = 1000
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
