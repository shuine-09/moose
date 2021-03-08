[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 900
    ny = 900
    xmax = 2.8
    ymax = 2.8
  []
[]

#[Adaptivity]
# initial_steps = 1
# max_h_level = 1
# stop_time = 1.0e-10
# initial_marker = err_bnds
#[./Markers]
#   [./err_bnds]
#     type = ErrorFractionMarker
#     coarsen = 0
#     refine = 0.9
#     indicator = ind_bnds
#   [../]
# [../]
# [./Indicators]
#    [./ind_bnds]
#      type = GradientJumpIndicator
#      variable = bnds
#   [../]
# [../]
#[]

[GlobalParams]
  displacements = 'u_x u_y'
[]

[UserObjects]
  [soln]
    type = SolutionUserObject
    mesh = b101_out.e
    timestep = 'LATEST'
    execute_on = 'initial'
  []
[]

[Variables]
  [d]
  []
[]

[AuxVariables]
  [global_strain]
    order = THIRD
    family = SCALAR
  []
  [u_x]
  []
  [u_y]
  []
  [bnds]
    [InitialCondition]
      type = FunctionIC
      function = f_bnds
      variable = bnds
    []
  []
  [c]
    order = FIRST
    family = LAGRANGE
    [InitialCondition]
      type = SmoothCircleIC
      x1 = 1.4
      y1 = 1.4
      radius = 0.5
      invalue = 1.0
      outvalue = 0.0
      int_width = 0.05
    []
  []
  [bounds_dummy]
    order = FIRST
    family = LAGRANGE
  []
[]

[Functions]
  [tfunc]
    type = ParsedFunction
    value = '0*t'
  []
  # [pressure]
  #   type = PiecewiseLinear
  #   data_file = 'pressure.csv'
  #   format = columns
  # []
  [pressure]
    type = PiecewiseLinear
    x = '0.0 10'
    y = '0.0 200'
  []
  # [./pressure]
  #   type = ParsedFunction
  #   value = '100'
  # [../]
  [f_bnds]
    type = SolutionFunction
    solution = soln
    from_variable = 'bnds'
  []
[]

[Bounds]
  [d_upper_bound]
    type = ConstantBoundsAux
    variable = bounds_dummy
    bounded_variable = d
    bound_type = upper
    bound_value = 1.0
  []
  [d_lower_bound]
    type = VariableOldValueBoundsAux
    variable = bounds_dummy
    bounded_variable = d
    bound_type = lower
  []
[]

[Kernels]
  [pfbulk]
    type = AllenCahn
    variable = d
    mob_name = L
    f_name = F
  []
  [dcdt]
    type = TimeDerivative
    variable = d
  []
  [acint]
    type = ACInterface
    variable = d
    mob_name = L
    kappa_name = kappa_op
  []
[]

[UserObjects]
  [global_strain_uo]
    type = GlobalStrainUserObject
    applied_stress_tensor = '-60 -60 -60 0 0 0'
    execute_on = 'Initial Linear Nonlinear'
  []
[]

[Materials]
  [pfbulkmat]
    type = GenericConstantMaterial
    prop_names = 'l visco'
    prop_values = '0.01 1e-3'
  []

  [pressure]
    type = GenericFunctionMaterial
    block = 0
    prop_names = fracture_pressure
    prop_values = pressure
  []

  [pressure_void]
    type = ParsedMaterial
    block = 0
    f_name = pressure_void
    args = 'c'
    material_property_names = 'pressure'
    function = 'pressure * c'
  []

  [gc]
    type = ParsedMaterial
    f_name = gc_prop
    #function = 'if(bnds < 0.75, if(bnds>0.25, 0.5, 2.5), 2.5)'
    function = 'if(bnds < 0.75 & c < 0.5, 0.0012, 0.012)'
    args = 'bnds c'
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
  [strain]
    type = ComputeSmallStrain
    global_strain = global_strain
  []
  [global_strain]
    type = ComputeGlobalStrain
    scalar_global_strain = global_strain
    global_strain_uo = global_strain_uo
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
      variable = 'd'
    []
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
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -ksp_type -snes_type -pc_factor_shift_type -pc_factor_shift_amount '
  petsc_options_value = 'lu preonly  vinewtonrsls NONZERO 1e-10'

  #  solve_type = PJFNK
  #  petsc_options_iname = '-pc_type -sub_pc_type -snes_type'
  #  petsc_options_value = 'asm lu vinewtonrsls'
  nl_rel_tol = 1e-6 ##nonlinear relative tolerance
  nl_abs_tol = 1e-6
  l_max_its = 10 ##max linear iterations Previous:200
  nl_max_its = 20 ##max nonlinear iterations Previous:50
  start_time = 0
  line_search = 'none'
  end_time = 2000
  dtmax = 1
  dtmin = 1e-14
  automatic_scaling = true
  #  [./TimeStepper]
  #    type = IterationAdaptiveDT
  #    dt = 1
  #    optimal_iterations = 10
  #    iteration_window = 0
  #    growth_factor = 1.2
  #    cutback_factor = 0.5
  #  [../]
[]

[Outputs]
  print_linear_converged_reason = false
  print_nonlinear_converged_reason = false
  print_linear_residuals = false

  #gnuplot = true
[]
