[GlobalParams]
  displacements = 'disp_x disp_y'
  #out_of_plane_strain = strain_zz
[]

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 16
    ny = 16
    xmin = 0
    xmax = 4
    ymin = 0
    ymax = 4
  []
  #uniform_refine = 4
[]

[Variables]
  # [./strain_zz]
  #   order = FIRST
  #   family = LAGRANGE
  # [../]
[]

[Adaptivity]
  steps = 1
  marker = box
  max_h_level = 8
  initial_steps = 8
  stop_time = 1.0e-10
  [./Markers]
    [./box]
      bottom_left = '1.6 1.8 0'
      inside = refine
      top_right = '2.4 2.2 0'
      outside = do_nothing
      type = BoxMarker
    [../]
  [../]
[]

[ICs]
  [./ini_c]
    type = PhaseFieldDamageIC
    variable = c
    l = 0.025
    d0 = 1.0
    x1 = 1.8
    x2 = 2.2
    y1 = 2
    y2 = 2
    z1 = 0
    z2 = 0
  [../]
[]

[Modules]
  [./TensorMechanics]
    [./Master]
      [./mech]
        add_variables = true
        strain = SMALL
        additional_generate_output = 'stress_yy stress_zz'
        #planar_formulation = WEAK_PLANE_STRESS
      [../]
    [../]
  [../]
  [./PhaseField]
    [./Nonconserved]
      [./c]
        free_energy = F
        kappa = kappa_op
        mobility = L
      [../]
    [../]
  [../]
[]

[AuxVariables]
  [./resid_x]
  [../]
  [./resid_y]
  [../]
  [./bounds_dummy]
    order = FIRST
    family = LAGRANGE
  [../]
  [./cod]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./elastic_energy]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./cod]
    type = PFCrackOpeningDisplacement
    displacements = 'disp_x disp_y'
    c = c
    variable = cod
  [../]
  [./elastic_energy]
    type = MaterialRealAux
    property = elastic_energy
    variable = elastic_energy
  [../]
[]

[BCs]
  [./yfix]
    type = DirichletBC
    preset = true
    variable = disp_y
    boundary = 'left right top bottom'
    value = 0
  [../]
  [./xfix]
    type = DirichletBC
    preset = true
    variable = disp_x
    boundary = 'left right top bottom'
    value = 0
  [../]
[]

[Bounds]
  [./c_upper_bound]
    type = ConstantBoundsAux
    variable = bounds_dummy
    bounded_variable = c
    bound_type = upper
    bound_value = 1.0
  [../]
  [./c_lower_bound]
    type = VariableOldValueBoundsAux
    variable = bounds_dummy
    bounded_variable = c
    bound_type = lower
  [../]
[]

[Functions]
  [./fracture_pressure]
    type = ParsedFunction
    value = '10*t'
  [../]
[]

[Materials]
  [./pressure]
    type = GenericFunctionMaterial
    prop_names = fracture_pressure
    prop_values = fracture_pressure
  [../]
  [./pfbulkmat]
    type = GenericConstantMaterial
    prop_names = 'gc_prop l visco'
    prop_values = '1 0.025 1e-6'
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
    function = 'gc_prop * l'
  [../]
  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1
    poissons_ratio = 0.2
  []
  [./damage_stress]
    type = ComputeLinearElasticPFFractureStress
    c = c
    E_name = 'elastic_energy'
    D_name = 'degradation'
    I_name = 'indicator_function'
    F_name = 'local_fracture_energy'
    decomposition_type = none
    use_current_history_variable = true
    use_snes_vi_solver = true
  [../]
  [./degradation]
    type = DerivativeParsedMaterial
    f_name = degradation
    args = 'c'
    function = '(1.0-c)^2*(1.0 - eta) + eta'
    constant_names       = 'eta'
    constant_expressions = '1e-12'
    derivative_order = 2
  [../]
  [./indicator_function]
    type = DerivativeParsedMaterial
    f_name = indicator_function
    args = 'c'
    function = 'c'
    derivative_order = 1
  [../]
  [./local_fracture_energy]
    type = DerivativeParsedMaterial
    f_name = local_fracture_energy
    args = 'c'
    material_property_names = 'gc_prop l'
    function = 'c^2 * gc_prop / 2 / l'
    derivative_order = 2
  [../]
  [./fracture_driving_energy]
    type = DerivativeSumMaterial
    args = c
    sum_materials = 'elastic_energy local_fracture_energy'
    derivative_order = 2
    f_name = F
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[VectorPostprocessors]
  [./cod1]
    type = LineValueSampler
    start_point = '1.5 4 0'
    end_point = '1.5 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod2]
    type = LineValueSampler
    start_point = '1.55 4 0'
    end_point = '1.55 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod3]
    type = LineValueSampler
    start_point = '1.6 4 0'
    end_point = '1.6 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod4]
    type = LineValueSampler
    start_point = '1.65 4 0'
    end_point = '1.65 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod5]
    type = LineValueSampler
    start_point = '1.7 4 0'
    end_point = '1.7 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod6]
    type = LineValueSampler
    start_point = '1.71 4 0'
    end_point = '1.71 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod7]
    type = LineValueSampler
    start_point = '1.72 4 0'
    end_point = '1.72 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod8]
    type = LineValueSampler
    start_point = '1.73 4 0'
    end_point = '1.73 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod9]
    type = LineValueSampler
    start_point = '1.74 4 0'
    end_point = '1.74 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod10]
    type = LineValueSampler
    start_point = '1.75 4 0'
    end_point = '1.75 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod11]
    type = LineValueSampler
    start_point = '1.76 4 0'
    end_point = '1.76 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod12]
    type = LineValueSampler
    start_point = '1.77 4 0'
    end_point = '1.77 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod13]
    type = LineValueSampler
    start_point = '1.78 4 0'
    end_point = '1.78 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod14]
    type = LineValueSampler
    start_point = '1.79 4 0'
    end_point = '1.79 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod15]
    type = LineValueSampler
    start_point = '1.8 4 0'
    end_point = '1.8 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod16]
    type = LineValueSampler
    start_point = '1.81 4 0'
    end_point = '1.81 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod17]
    type = LineValueSampler
    start_point = '1.82 4 0'
    end_point = '1.82 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod18]
    type = LineValueSampler
    start_point = '1.83 4 0'
    end_point = '1.83 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod19]
    type = LineValueSampler
    start_point = '1.84 4 0'
    end_point = '1.84 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod20]
    type = LineValueSampler
    start_point = '1.85 4 0'
    end_point = '1.85 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod21]
    type = LineValueSampler
    start_point = '1.86 4 0'
    end_point = '1.86 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod22]
    type = LineValueSampler
    start_point = '1.87 4 0'
    end_point = '1.87 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod23]
    type = LineValueSampler
    start_point = '1.88 4 0'
    end_point = '1.88 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod24]
    type = LineValueSampler
    start_point = '1.89 4 0'
    end_point = '1.89 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod25]
    type = LineValueSampler
    start_point = '1.9 4 0'
    end_point = '1.9 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod26]
    type = LineValueSampler
    start_point = '1.91 4 0'
    end_point = '1.91 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod27]
    type = LineValueSampler
    start_point = '1.92 4 0'
    end_point = '1.92 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod28]
    type = LineValueSampler
    start_point = '1.93 4 0'
    end_point = '1.93 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod29]
    type = LineValueSampler
    start_point = '1.94 4 0'
    end_point = '1.94 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod30]
    type = LineValueSampler
    start_point = '1.95 4 0'
    end_point = '1.95 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod31]
    type = LineValueSampler
    start_point = '1.96 4 0'
    end_point = '1.96 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod32]
    type = LineValueSampler
    start_point = '1.97 4 0'
    end_point = '1.97 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod33]
    type = LineValueSampler
    start_point = '1.98 4 0'
    end_point = '1.98 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod34]
    type = LineValueSampler
    start_point = '1.99 4 0'
    end_point = '1.99 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod35]
    type = LineValueSampler
    start_point = '2.0 4 0'
    end_point = '2.0 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod36]
    type = LineValueSampler
    start_point = '2.01 4 0'
    end_point = '2.01 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod37]
    type = LineValueSampler
    start_point = '2.02 4 0'
    end_point = '2.02 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod38]
    type = LineValueSampler
    start_point = '2.03 4 0'
    end_point = '2.03 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod39]
    type = LineValueSampler
    start_point = '2.04 4 0'
    end_point = '2.04 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod40]
    type = LineValueSampler
    start_point = '2.05 4 0'
    end_point = '2.05 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod41]
    type = LineValueSampler
    start_point = '2.06 4 0'
    end_point = '2.06 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod42]
    type = LineValueSampler
    start_point = '2.07 4 0'
    end_point = '2.07 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod43]
    type = LineValueSampler
    start_point = '2.08 4 0'
    end_point = '2.08 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod44]
    type = LineValueSampler
    start_point = '2.09 4 0'
    end_point = '2.09 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod45]
    type = LineValueSampler
    start_point = '2.1 4 0'
    end_point = '2.1 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod46]
    type = LineValueSampler
    start_point = '2.11 4 0'
    end_point = '2.11 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod47]
    type = LineValueSampler
    start_point = '2.12 4 0'
    end_point = '2.12 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod48]
    type = LineValueSampler
    start_point = '2.13 4 0'
    end_point = '2.13 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod49]
    type = LineValueSampler
    start_point = '2.14 4 0'
    end_point = '2.14 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod50]
    type = LineValueSampler
    start_point = '2.15 4 0'
    end_point = '2.15 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod51]
    type = LineValueSampler
    start_point = '2.16 4 0'
    end_point = '2.16 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod52]
    type = LineValueSampler
    start_point = '2.17 4 0'
    end_point = '2.17 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod53]
    type = LineValueSampler
    start_point = '2.18 4 0'
    end_point = '2.18 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod54]
    type = LineValueSampler
    start_point = '2.19 4 0'
    end_point = '2.19 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod55]
    type = LineValueSampler
    start_point = '2.20 4 0'
    end_point = '2.20 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod56]
    type = LineValueSampler
    start_point = '2.21 4 0'
    end_point = '2.21 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod57]
    type = LineValueSampler
    start_point = '2.22 4 0'
    end_point = '2.22 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod58]
    type = LineValueSampler
    start_point = '2.23 4 0'
    end_point = '2.23 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod59]
    type = LineValueSampler
    start_point = '2.24 4 0'
    end_point = '2.24 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod60]
    type = LineValueSampler
    start_point = '2.25 4 0'
    end_point = '2.25 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod61]
    type = LineValueSampler
    start_point = '2.26 4 0'
    end_point = '2.26 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod62]
    type = LineValueSampler
    start_point = '2.27 4 0'
    end_point = '2.27 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod63]
    type = LineValueSampler
    start_point = '2.28 4 0'
    end_point = '2.28 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod64]
    type = LineValueSampler
    start_point = '2.29 4 0'
    end_point = '2.29 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod65]
    type = LineValueSampler
    start_point = '2.3 4 0'
    end_point = '2.3 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod66]
    type = LineValueSampler
    start_point = '2.35 4 0'
    end_point = '2.35 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod67]
    type = LineValueSampler
    start_point = '2.4 4 0'
    end_point = '2.4 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod68]
    type = LineValueSampler
    start_point = '2.45 4 0'
    end_point = '2.45 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
  [./cod69]
    type = LineValueSampler
    start_point = '2.5 4 0'
    end_point = '2.5 0 0'
    variable = 'cod'
    num_points = 10000
    sort_by = 'y'
  [../]
[]

[Executioner]
  type = Transient

  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -snes_type'
  petsc_options_value = 'lu       superlu_dist vinewtonrsls'

  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-8
  l_max_its = 10
  nl_max_its = 1000
  line_search = 'none'

  dt = 1e-4
  end_time = 1e-4
  #automatic_scaling = true
[]

[Outputs]
  exodus = true
  csv = true
[]
