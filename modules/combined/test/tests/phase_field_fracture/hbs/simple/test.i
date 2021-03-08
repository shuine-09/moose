[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 200
    ny = 200
    xmax = 2.8
    ymax = 2.8
  []
  [cnode]
    type = ExtraNodesetGenerator
    coord = '0 0 0'
    new_boundary = 100
    input = gen
  []
[]

[GlobalParams]
  displacements = 'u_x u_y'
[]

[Variables]
  [global_strain]
    order = THIRD
    family = SCALAR
  []
[]

[Modules]
  [TensorMechanics]
    [Master]
      [mech]
        add_variables = true
        strain = SMALL
        additional_generate_output = 'stress_yy stress_xy stress_xx strain_xx strain_xy strain_yy '
                                     'strain_zz hydrostatic_stress mid_principal_stress '
                                     'min_principal_stress max_principal_stress'
        decomposition_method = EigenSolution
        global_strain = global_strain
        save_in = 'force_x force_y'
      []
    []
    [GlobalStrain]
      [global_strain]
        scalar_global_strain = global_strain
        displacements = 'u_x u_y'
        auxiliary_displacements = 'disp_x disp_y'
        global_displacements = 'ug_x ug_y'
      []
    []
  []
[]

[AuxVariables]
  [force_x]
    order = FIRST
    family = LAGRANGE
  []
  [force_y]
    order = FIRST
    family = LAGRANGE
  []
  [c]
    order = FIRST
    family = LAGRANGE
    [InitialCondition]
      type = SpecifiedSmoothCircleIC
      invalue = 1.0
      outvalue = 0
      int_width = 0.25
      x_positions = '1.4'
      z_positions = '0 0'
      y_positions = '1.4'
      radii = '0.25'
    []
  []
  [bounds_dummy]
    order = FIRST
    family = LAGRANGE
  []
[]

[Functions]
  [pressure]
    type = PiecewiseLinear
    x = '0.0 10'
    y = '0.0 10'
  []
[]

[Materials]
  [pressure]
    type = GenericFunctionMaterial
    block = 0
    prop_names = fracture_pressure
    prop_values = pressure
  []
  [stress_void]
    type = ComputeLinearElasticStress
    block = 0
    base_name = void
  []
  [strain_void]
    type = ComputeSmallStrain
    block = 0
    base_name = void
  []

  [strain]
    type = ComputeSmallStrain
    block = 0
    base_name = matrix
  []

  [damage_stress]
    type = ComputeLinearElasticStress
    base_name = matrix
  []
  [const_stress]
    type = ComputeExtraStressConstant
    block = 0
    base_name = void
    extra_stress_tensor = '-1 -1 -1 0 0 0'
    prefactor = fracture_pressure
  []
  [global_stress]
    type = TwoPhaseStressMaterial
    base_A = matrix
    base_B = void
  []
  [switching]
    type = SwitchingFunctionMaterial
    eta = c
  []
  [elasticity_tensor]
    type = ComputeElasticityTensor
    C_ijkl = '385000 0.23'
    fill_method = symmetric_isotropic_E_nu
    base_name = matrix
  []
  [elasticity_tensor_void]
    type = ComputeElasticityTensor
    C_ijkl = '3.85 0.23'
    fill_method = symmetric_isotropic_E_nu
    base_name = void
  []
[]

[BCs]
  [Periodic]
    [all]
      auto_direction = 'x y'
      variable = 'u_x u_y'
    []
  []
  # [./yfix]
  #   type = DirichletBC
  #   variable = u_y
  #   boundary = 100
  #   value = 0
  # [../]
  # [./xfix]
  #   type = DirichletBC
  #   variable = u_x
  #   boundary = 100
  #   value = 0
  # [../]
  # [./Pressure]
  #   [./Pressure]
  #     boundary = 'top right'
  #     factor = 30
  #     function = 1
  #   [../]
  # [../]
[]

[Postprocessors]
  [ave_stress_top]
    type = SideAverageValue
    variable = stress_yy
    boundary = top
  []
  [disp_y_top]
    type = SideAverageValue
    variable = disp_y
    boundary = top
  []
  [react_y_top]
    type = NodalSum
    variable = force_y
    boundary = top
  []
[]

[Preconditioning]
  active = 'smp'
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  # petsc_options_iname = '-pc_type -sub_pc_type -snes_type'
  # petsc_options_value = 'asm lu vinewtonrsls'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-6 ##nonlinear relative tolerance
  nl_abs_tol = 1e-6
  l_max_its = 10 ##max linear iterations Previous:200
  nl_max_its = 20 ##max nonlinear iterations Previous:50
  start_time = 0
  line_search = 'none'
  end_time = 2000
  num_steps = 10
  dtmax = 1
  dtmin = 1e-14
  automatic_scaling = true
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1
    optimal_iterations = 10
    iteration_window = 0
    growth_factor = 1.2
    cutback_factor = 0.5
  []

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
