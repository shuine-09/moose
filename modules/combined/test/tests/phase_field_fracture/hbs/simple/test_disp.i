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
    coord = '1.4 1.4 0'
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

[Adaptivity]
  marker = EFM_0
  max_h_level = 1
  #cycles_per_step = 1
  [./Markers]
    [./EFM_0]
      type = ErrorFractionMarker
      coarsen = 0.05
      refine = 0.95
      indicator = GJI_0
    [../]
  [../]
  [./Indicators]
    [./GJI_0]
      type = GradientJumpIndicator
      variable = d
    [../]
  [../]
[]


[MultiApps]
  [damage]
    type = TransientMultiApp
    input_files = 'test_c.i'
  []
[]

[Transfers]
  [to_disp_x]
    type = MultiAppCopyTransfer
    multi_app = 'damage'
    direction = to_multiapp
    source_variable = 'u_x'
    variable = 'u_x'
  []
  [to_disp_y]
    type = MultiAppCopyTransfer
    multi_app = 'damage'
    direction = to_multiapp
    source_variable = 'u_y'
    variable = 'u_y'
  []
  [to_global_strain]
    type = MultiAppScalarToAuxScalarTransfer
    multi_app = 'damage'
    direction = to_multiapp
    source_variable = 'global_strain'
    to_aux_scalar = 'global_strain'
  []
  [from_c]
    type = MultiAppCopyTransfer
    multi_app = 'damage'
    direction = from_multiapp
    source_variable = 'd'
    variable  = 'd'
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
  [d]
  []
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
      int_width = 0.05
      x_positions = '1.4'
      z_positions = '0 0'
      y_positions = '1.4'
      radii = '0.5'
    []
  []
  # [gc]
  #   order = FIRST
  #   family = LAGRANGE
  #   [InitialCondition]
  #     type = VolumeWeightedWeibull
  #     reference_volume = 0.000196 #This is the volume of an element for a 100x100 mesh
  #     weibull_modulus = 15.0
  #     median = 0.0012
  #   []
  # []
[]

[Functions]
  [pressure]
    type = PiecewiseLinear
    x = '0.0 10'
    y = '0.0 200'
  []
[]

[Functions]
  [./func]
    type = ParsedFunction
    value = 'if((y<1.42 & y>1.38) | (x<1.42 & x>1.38), 0.0012, 1.2)'
  [../]
[]

[Materials]
  [pressure]
    type = GenericFunctionMaterial
    block = 0
    prop_names = pressure
    prop_values = pressure
  []
  [pressure_void]
    type = ParsedMaterial
    block = 0
    f_name = pressure_void
    args = 'c'
    material_property_names = 'pressure'
    function = 'pressure * c'
    outputs = exodus
  []

  # [gc_prop]
  #   type = ParsedMaterial
  #   block = 0
  #   f_name = gc_prop
  #   args = 'gc'
  #   function = 'gc'
  #   outputs = exodus
  # []

  [gc_prop]
    type = GenericFunctionMaterial
    prop_names = gc_prop
    prop_values = func
    outputs = exodus
  []

  [pfbulkmat]
    type = GenericConstantMaterial
    prop_names = 'l visco'
    prop_values = '0.021 1e-3'
  []
  [define_mobility]
    type = ParsedMaterial
    material_property_names = 'gc_prop visco'
    f_name = L
    function = '1.0/(gc_prop * visco)'
  []
  [define_kappa]
    type = ParsedMaterial
    material_property_names = 'gc_prop l'
    f_name = kappa_op
    function = 'gc_prop * l / 3.14159 * 2'
  []

  [damage_stress]
    type = ComputeLinearElasticPFFractureStress
    c = d
    E_name = 'elastic_energy'
    D_name = 'degradation'
    F_name = 'local_fracture_energy'
    I_name = 'indicator_function'
    decomposition_type = none
    use_snes_vi_solver = true
    output_properties = 'hist'
    outputs = exodus
  []
  [indicator_function]
    type = DerivativeParsedMaterial
    f_name = indicator_function
    args = 'd'
    function = 'd*d'
    derivative_order = 2
  []
  [degradation]
    type = DerivativeParsedMaterial
    f_name = degradation
    args = 'd'
    function = '((1.0-d)^2+eta)/((1.0-d)^2+d*(1-0.5*d)*(4/3.14159/l*E*gc_prop/sigma^2))'
    material_property_names = 'gc_prop l'
    constant_names = 'E sigma eta'
    constant_expressions = '385000 130 1e-4'
    derivative_order = 2
  []
  [fracture_energy]
    type = DerivativeParsedMaterial
    f_name = local_fracture_energy
    args = 'd'
    material_property_names = 'gc_prop l'
    function = 'gc_prop/l/3.14159*(2*d-d^2)'
    derivative_order = 2
  []
  [fracture_driving_energy]
    type = DerivativeParsedMaterial
    args = 'c d'
    material_property_names = 'elastic_energy(d) local_fracture_energy(d)'
    function = '(1-c)*elastic_energy + local_fracture_energy'
    derivative_order = 2
    f_name = F
  []

  [const_stress]
    type = ComputeExtraStressConstant
    block = 0
    extra_stress_tensor = '-1 -1 -1 0 0 0'
    prefactor = pressure_void
  []

  # [elasticity_tensor]
  #   type = ComputeElasticityTensor
  #   block = 0
  #   C_ijkl = '385000 0.23'
  #   fill_method = symmetric_isotropic_E_nu
  # []
  [elasticity_tensor]
    type = ComputeConcentrationDependentElasticityTensor
    block = 0
    c = c
    C1_ijkl = '3.85 0.23'
    C0_ijkl = '385000 0.23'
    fill_method1 = symmetric_isotropic_E_nu
    fill_method0 = symmetric_isotropic_E_nu
  []
[]

[BCs]
  [Periodic]
    [all]
      auto_direction = 'x y'
      variable = 'u_x u_y'
    []
  []
  [yfix]
    type = DirichletBC
    variable = u_y
    boundary = 100
    value = 0
  []
  [xfix]
    type = DirichletBC
    variable = u_x
    boundary = 100
    value = 0
  []
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
  nl_rel_tol = 1e-10 ##nonlinear relative tolerance
  nl_abs_tol = 1e-10
  l_max_its = 10 ##max linear iterations Previous:200
  nl_max_its = 20 ##max nonlinear iterations Previous:50
  start_time = 0
  line_search = 'none'
  end_time = 2000
  num_steps = 100
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

  picard_max_its = 20
  picard_rel_tol = 1e-10
  picard_abs_tol = 1e-10
  accept_on_max_picard_iteration = true
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
