[GlobalParams]
  displacements = 'disp_x disp_y'
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
  [./c]
  [../]
[]

[Adaptivity]
  steps = 1
  marker = box
  max_h_level = 5
  initial_steps = 5
  stop_time = 1.0e-10
  [./Markers]
    [./box]
      bottom_left = '1.2 1.8 0'
      inside = refine
      top_right = '2.8 2.2 0'
      outside = do_nothing
      type = BoxMarker
    [../]
  [../]
[]

[ICs]
  # [./ini_c]
  #   type = PhaseFieldDamageIC
  #   variable = c
  #   l = 0.06
  #   d0 = 1.0
  #   x1 = 1.8
  #   x2 = 2.2
  #   y1 = 2
  #   y2 = 2
  #   z1 = 0
  #   z2 = 0
  # [../]
  [./ini_c]
    type = BoundingBoxIC
    variable = c
    x1 = 1.8
    y1 = 1.99999
    x2 = 2.2
    y2 = 2.00001
    inside = 1.0
    outside = 0.0
  [../]
[]

[Modules]
  [./TensorMechanics]
    [./Master]
      [./mech]
        add_variables = true
        strain = SMALL
        additional_generate_output = 'stress_yy stress_zz'
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

[Kernels]
  [./pf_c]
    type = PhaseFieldPressurizedFractureDamage
    variable = c
    displacements = 'disp_x disp_y'
  [../]
  [./pf_disp_x]
    type = PhaseFieldPressurizedFractureMechanics
    variable = disp_x
    component = 0
    c = c
  [../]
  [./pf_disp_y]
    type = PhaseFieldPressurizedFractureMechanics
    variable = disp_y
    component = 1
    c = c
  [../]
  [./ACBulk]
    type = AllenCahn
    variable = c
    f_name = F
  [../]

  [./ACInterface]
    type = ACInterface
    variable = c
    kappa_name = kappa_op
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
    value = 't'
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
    prop_values = '0.001 0.04 1e-6'
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
    F_name = 'local_fracture_energy'
    I_name = 'indicator_function'
    decomposition_type = none
    use_current_history_variable = true
    use_snes_vi_solver = true
  [../]
  [./indicator_function]
    type = DerivativeParsedMaterial
    f_name = indicator_function
    args = 'c'
    function = 'c*c'
    derivative_order = 1
  [../]
  [./degradation]
    type = DerivativeParsedMaterial
    f_name = degradation
    args = 'c'
    function = '(1.0-c)^2*(1.0 - eta) + eta'
    constant_names       = 'eta'
    constant_expressions = '1e-8'
    derivative_order = 2
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

[Postprocessors/crack_length]
  type = PFCrackVolume
  l = 0.04
  variable = c
[]

[Executioner]
  type = Transient

  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -snes_type'
  petsc_options_value = 'lu       superlu_dist vinewtonrsls'

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6
  l_max_its = 10
  nl_max_its = 50
  line_search = 'none'

  dt = 2.5e-4
  dtmin = 1e-7
  end_time = 350e-4
  #automatic_scaling = true
[]

[Outputs]
  file_base = sneddon_a02
  exodus = true
  csv = true
[]
