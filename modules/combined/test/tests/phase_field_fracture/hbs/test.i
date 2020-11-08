# [Mesh]
#   type = FileMesh
#   file = hbs_bub_test.e
# []

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmax = 950
    ymax = 950
    nx = 300
    ny = 300
  []
[]

[GlobalParams]
  op_num = 15
  var_name_base = gr
  displacements = 'disp_x disp_y'
[]

[UserObjects]
  [./soln]
    type = SolutionUserObject
    mesh = hbs_bub_test.e
    timestep = 'LATEST'
    execute_on = 'initial'
  [../]
[]

[Variables]
  [./d]
  [../]
  [./disp_x]
  [../]
  [./disp_y]
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
        save_in = 'force_x force_y'
      [../]
    [../]
  [../]
[]

[AuxVariables]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./force_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./force_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./bnds]
    [./InitialCondition]
      type = FunctionIC
      function = f_bnds
      variable = bnds
    [../]
  [../]
  [./c]
    [./InitialCondition]
      type = FunctionIC
      function = f_bubs
      variable = c
    [../]
  [../]
  [./gr0]
    [./InitialCondition]
      type = FunctionIC
      variable = gr0
      function = f_gr0
    [../]
  [../]
  [./gr1]
    [./InitialCondition]
      type = FunctionIC
      variable = gr1
      function = f_gr1
    [../]
  [../]
  [./gr2]
    [./InitialCondition]
      type = FunctionIC
      variable = gr2
      function = f_gr2
    [../]
  [../]
  [./gr3]
    [./InitialCondition]
      type = FunctionIC
      variable = gr3
      function = f_gr3
    [../]
  [../]
  [./gr4]
    [./InitialCondition]
      type = FunctionIC
      variable = gr4
      function = f_gr4
    [../]
  [../]
  [./gr5]
    [./InitialCondition]
      type = FunctionIC
      variable = gr5
      function = f_gr5
    [../]
  [../]
  [./gr6]
    [./InitialCondition]
      type = FunctionIC
      variable = gr6
      function = f_gr6
    [../]
  [../]
  [./gr7]
    [./InitialCondition]
      type = FunctionIC
      variable = gr7
      function = f_gr7
    [../]
  [../]
  [./gr8]
    [./InitialCondition]
      type = FunctionIC
      variable = gr8
      function = f_gr8
    [../]
  [../]
  [./gr9]
    [./InitialCondition]
      type = FunctionIC
      variable = gr9
      function = f_gr9
    [../]
  [../]
  [./gr10]
    [./InitialCondition]
      type = FunctionIC
      variable = gr10
      function = f_gr10
    [../]
  [../]
  [./gr11]
    [./InitialCondition]
      type = FunctionIC
      variable = gr11
      function = f_gr11
    [../]
  [../]
  [./gr12]
    [./InitialCondition]
      type = FunctionIC
      variable = gr12
      function = f_gr12
    [../]
  [../]
  [./gr13]
    [./InitialCondition]
      type = FunctionIC
      variable = gr13
      function = f_gr13
    [../]
  [../]
  [./gr14]
    [./InitialCondition]
      type = FunctionIC
      variable = gr14
      function = f_gr14
    [../]
  [../]

  # [./euler_angle]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  [./C1111]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  # [./var_indices]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  [./strain_yy]
    family = MONOMIAL
    order = CONSTANT
  [../]
  # [./unique_grains]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  [./bounds_dummy]
    order = FIRST
    family = LAGRANGE
  [../]
[]


[Functions]
  [./tfunc]
    type = ParsedFunction
    value = '0*t'
  [../]
 [./pressure]
   type = PiecewiseLinear
   data_file = 'bubble_pressure1.csv'
   format = columns
 [../]
  # [./pressure]
  #   type = ParsedFunction
  #   value = '100'
  # [../]
  [./f_bnds]
    type = SolutionFunction
    solution = soln
    from_variable = 'bnds'
  [../]
  [./f_bubs]
    type = SolutionFunction
    solution = soln
    from_variable = 'etab'
  [../]
  [./f_gr0]
    type = SolutionFunction
    solution = soln
    from_variable = 'eta0'
  [../]
  [./f_gr1]
    type = SolutionFunction
    solution = soln
    from_variable = 'eta1'
  [../]
  [./f_gr2]
    type = SolutionFunction
    solution = soln
    from_variable = 'eta2'
  [../]
  [./f_gr3]
    type = SolutionFunction
    solution = soln
    from_variable = 'eta3'
  [../]
  [./f_gr4]
    type = SolutionFunction
    solution = soln
    from_variable = 'eta4'
  [../]
  [./f_gr5]
    type = SolutionFunction
    solution = soln
    from_variable = 'eta5'
  [../]
  [./f_gr6]
    type = SolutionFunction
    solution = soln
    from_variable = 'eta6'
  [../]
  [./f_gr7]
    type = SolutionFunction
    solution = soln
    from_variable = 'eta7'
  [../]
  [./f_gr8]
    type = SolutionFunction
    solution = soln
    from_variable = 'eta8'
  [../]
  [./f_gr9]
    type = SolutionFunction
    solution = soln
    from_variable = 'eta9'
  [../]
  [./f_gr10]
    type = SolutionFunction
    solution = soln
    from_variable = 'eta10'
  [../]
  [./f_gr11]
    type = SolutionFunction
    solution = soln
    from_variable = 'eta11'
  [../]
  [./f_gr12]
    type = SolutionFunction
    solution = soln
    from_variable = 'eta12'
  [../]
  [./f_gr13]
    type = SolutionFunction
    solution = soln
    from_variable = 'eta13'
  [../]
  [./f_gr14]
    type = SolutionFunction
    solution = soln
    from_variable = 'eta14'
  [../]

[]

[Bounds]
  [./d_upper_bound]
    type = ConstantBoundsAux
    variable = bounds_dummy
    bounded_variable = d
    bound_type = upper
    bound_value = 1.0
  [../]
  [./d_lower_bound]
    type = VariableOldValueBoundsAux
    variable = bounds_dummy
    bounded_variable = d
    bound_type = lower
  [../]
[]

[Kernels]
  [./pfbulk]
    type = AllenCahn
    variable = d
    mob_name = L
    f_name = F
  [../]
  [./solid_x]
    type = PhaseFieldFractureMechanicsOffDiag
    variable = disp_x
    component = 0
    c = d
  [../]
  [./solid_y]
    type = PhaseFieldFractureMechanicsOffDiag
    variable = disp_y
    component = 1
    c = d
  [../]
  [./dcdt]
    type = TimeDerivative
    variable = d
  [../]
  [./acint]
    type = ACInterface
    variable = d
    mob_name = L
    kappa_name = kappa_op
  [../]
[]

[AuxKernels]
  [./bnds_aux]
    type = BndsCalcAux
    variable = bnds
    execute_on = timestep_end
  [../]
  [./C1111]
    type = RankFourAux
    variable = C1111
    rank_four_tensor = matrix_elasticity_tensor
    index_l = 0
    index_j = 0
    index_k = 0
    index_i = 0
    execute_on = timestep_end
  [../]
[]

[Materials]
  [./pfbulkmat]
    type = GenericConstantMaterial
    prop_names = 'l visco'
    prop_values = '10 1e-3'
  [../]
  [pressure]
    type = GenericFunctionMaterial
    block = 0
    prop_names = fracture_pressure
    prop_values = pressure
  []
  [./gc]
    type = ParsedMaterial
    f_name = gc_prop
    #function = 'if(bnds < 0.75, if(bnds>0.25, 0.5, 2.5), 2.5)'
    function = 'if(bnds < 0.75 & c < 0.5, 0.0012, 0.012)'
    args = 'bnds c'
  [../]

  [./define_mobility]
    type = ParsedMaterial
    material_property_names = 'gc_prop visco'
    f_name = L
    function = '1.0/(gc_prop * visco)'
  [../]
  [./define_kappa]
    type = ParsedMaterial
    material_property_names = 'gc_prop l'
    f_name = kappa_op
    function = 'gc_prop * l / 3.14159 * 2'
  [../]

  [./stress_void]
    type = ComputeLinearElasticStress
    block = 0
    base_name = void
  [../]
  [./strain_void]
    type = ComputeSmallStrain
    block = 0
    base_name = void
  [../]

  [./strain]
    type = ComputeSmallStrain
    block = 0
    base_name = matrix
  [../]

  [./damage_stress]
    type = ComputeLinearElasticPFFractureStress
    c = d
    E_name = 'elastic_energy'
    D_name = 'degradation'
    F_name = 'local_fracture_energy'
    I_name = 'indicator_function'
    decomposition_type = strain_vol_dev
    #decomposition_type = none
    use_snes_vi_solver = true
    base_name = matrix
  [../]
  [./indicator_function]
    type = DerivativeParsedMaterial
    f_name = indicator_function
    args = 'd'
    function = 'd*d'
    derivative_order = 2
  [../]
  [./degradation]
    type = DerivativeParsedMaterial
    f_name = degradation
    args = 'd'
    function = '((1.0-d)^2+eta)/((1.0-d)^2+d*(1-0.5*d)*(4/3.14159/l*E*gc_prop/sigma^2))'
    material_property_names = 'gc_prop l'
    constant_names       = 'E sigma eta'
    constant_expressions = '385000 130 1e-4'
    derivative_order = 2
  [../]
  [./fracture_energy]
    type = DerivativeParsedMaterial
    f_name = local_fracture_energy
    args = 'd'
    material_property_names = 'gc_prop l'
    function = 'gc_prop/l/3.14159*(2*d-d^2)'
    derivative_order = 2
  [../]
  [./fracture_driving_energy]
    type = DerivativeSumMaterial
    args = d
    sum_materials = 'elastic_energy local_fracture_energy'
    derivative_order = 2
    f_name = F
  [../]
  [./const_stress]
    type = ComputeExtraStressConstant
    block = 0
    base_name = void
    extra_stress_tensor = '-1 -1 -1 0 0 0'
    prefactor = fracture_pressure
  [../]
  [./global_stress]
    type = TwoPhaseStressMaterial
    base_A = matrix
    base_B = void
  [../]
  [./switching]
    type = SwitchingFunctionMaterial
    eta = c
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    C_ijkl = '385000 0.23'
    fill_method = symmetric_isotropic_E_nu
    base_name = matrix
  [../]
  [./elasticity_tensor_void]
    type = ComputeElasticityTensor
    C_ijkl = '3.85 0.23'
    fill_method = symmetric_isotropic_E_nu
    base_name = void
  [../]
[]

[BCs]
  [./yfix]
    type = PresetBC
    variable = disp_y
    boundary = 'bottom'
    value = 0
  [../]
  [./xfix]
    type = PresetBC
    variable = disp_x
    boundary = 'left'
    value = 0
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
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type -snes_type'
  petsc_options_value = 'lu superlu_dist vinewtonrsls'
  nl_rel_tol = 1e-6  ##nonlinear relative tolerance
  nl_abs_tol = 1e-6
  l_max_its = 10   ##max linear iterations Previous:200
  nl_max_its = 15  ##max nonlinear iterations Previous:50
  start_time=0
  line_search = 'none'
  end_time = 2000
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
#num_grids = 1
[]

[Outputs]
  [exodus]
    type = Exodus
    interval = 1
  []
  csv = true
#gnuplot = true
[]
