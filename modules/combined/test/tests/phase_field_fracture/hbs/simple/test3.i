[Mesh]
  [file]
    type = FileMeshGenerator
    file = hole.e
  []
[]

[GlobalParams]
  displacements = 'u_x u_y'
[]

[Variables]
  [./global_strain]
    order = THIRD
    family = SCALAR
  [../]
[]

[Modules]
  [./TensorMechanics]
    [./Master]
      [./mech]
        add_variables = true
        strain = SMALL
        additional_generate_output = 'stress_yy stress_xy stress_xx strain_xx strain_xy strain_yy strain_zz hydrostatic_stress mid_principal_stress min_principal_stress max_principal_stress'
        decomposition_method = EigenSolution
        global_strain = global_strain
        save_in = 'force_x force_y'
      [../]
    [../]
    [./GlobalStrain]
      [./global_strain]
        scalar_global_strain = global_strain
        displacements = 'u_x u_y'
        auxiliary_displacements = 'disp_x disp_y'
        global_displacements = 'ug_x ug_y'
      [../]
    [../]
  [../]
[]

[AuxVariables]
  [./force_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./force_y]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Materials]
  [./stress]
    type = ComputeLinearElasticStress
    block = 0
  [../]

  [./elasticity_tensor]
    type = ComputeElasticityTensor
    C_ijkl = '385000 0.23'
    fill_method = symmetric_isotropic_E_nu
  [../]
[]

[BCs]
  [./Periodic]
    [./all]
      auto_direction = 'x y'
      variable = 'u_x u_y'
    [../]
  [../]
  [./yfix]
    type = DirichletBC
    variable = u_y
    boundary = 100
    value = 0
  [../]
  [./xfix]
    type = DirichletBC
    variable = u_x
    boundary = 100
    value = 0
  [../]
  [./Pressure]
    [./Pressure]
      boundary = 'top right'
      factor = 30
      function = 1
    [../]
  [../]
[]

[Functions]
  [./pressure]
    type = PiecewiseLinear
    x = '0 1.0 2.0'
    y = '0 100 100'
  [../]
[]


[Postprocessors]
  [./ave_stress_top]
    type = SideAverageValue
    variable = stress_yy
    boundary = top
  [../]
  [./disp_y_top]
    type = SideAverageValue
    variable = disp_y
    boundary = top
  [../]
  [./react_y_top]
    type = NodalSum
    variable = force_y
    boundary = top
  [../]
[]

[Preconditioning]
  active = 'smp'
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  # petsc_options_iname = '-pc_type -sub_pc_type -snes_type'
  # petsc_options_value = 'asm lu vinewtonrsls'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-6  ##nonlinear relative tolerance
  nl_abs_tol = 1e-6
  l_max_its = 10   ##max linear iterations Previous:200
  nl_max_its = 20  ##max nonlinear iterations Previous:50
  start_time=0
  line_search = 'none'
  end_time = 2000
  num_steps = 10
  dtmax = 1
  dtmin = 1e-14
  automatic_scaling = true
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1
    optimal_iterations = 10
    iteration_window = 0
    growth_factor = 1.2
    cutback_factor = 0.5
  [../]

  # picard_max_its = 20
  # picard_rel_tol = 1e-6
  # picard_abs_tol = 1e-6
  # accept_on_max_picard_iteration = true
[]

[Outputs]
  [exodus]
    type = Exodus
    interval = 1
    execute_on = 'initial timestep_end'
  []
  csv = true
#gnuplot = true
[]
