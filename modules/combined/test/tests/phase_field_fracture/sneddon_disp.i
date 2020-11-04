a = 0.3
end = '${fparse 1.2*sqrt(1.0*0.001/(1-0.2^2)/pi/a)}'
W = '${fparse 4/0.2*a}'
H = '${fparse 4/0.2*a}'
H2 = '${fparse 2/0.2*a}'
nx = '${fparse ceil(16/0.2*a)}'
ny = '${fparse ceil(8/0.2*a)}'
a_left = '${fparse W/2.0-1}'
a_right = '${fparse W/2.0+1}'
a_top = '${fparse H2+0.2}'
a_bottom = '${fparse H2-0.2}'
x1 = '${fparse H2-a}'
x2 = '${fparse H2+a}'
y1 = '${fparse H2-0.0001}'
y2 = '${fparse H2+0.0001}'

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[MultiApps]
  [damage]
    type = TransientMultiApp
    input_files = 'sneddon_c.i'
  []
[]

[Transfers]
  [to_disp_x]
    type = MultiAppCopyTransfer
    multi_app = 'damage'
    direction = to_multiapp
    source_variable = 'disp_x'
    variable = 'disp_x'
  []
  [to_disp_y]
    type = MultiAppCopyTransfer
    multi_app = 'damage'
    direction = to_multiapp
    source_variable = 'disp_y'
    variable = 'disp_y'
  []
  [from_c]
    type = MultiAppCopyTransfer
    multi_app = 'damage'
    direction = from_multiapp
    source_variable = 'c'
    variable = 'c'
  []
[]

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmax = ${W}
    ymax = ${H}
    ymin = ${H2}
    nx = ${nx}
    ny = ${ny}
  []
[]

[AuxVariables]
  [c]
  []
[]

[Adaptivity]
  marker = box
  max_h_level = 5
  initial_steps = 5
  [Markers]
    [box]
      type = BoxMarker
      bottom_left = '${a_left} ${a_bottom} 0'
      top_right = '${a_right} ${a_top} 0'
      inside = refine
      outside = do_nothing
    []
  []
[]

[ICs]
  [ini_c]
    type = BoundingBoxIC
    variable = c
    x1 = ${x1}
    y1 = ${y1}
    x2 = ${x2}
    y2 = ${y2}
    inside = 1.0
    outside = 0.0
  []
[]

[Modules]
  [TensorMechanics]
    [Master]
      [mech]
        add_variables = true
        strain = SMALL
        additional_generate_output = 'stress_yy stress_zz'
      []
    []
  []
[]

[BCs]
  [yfix]
    type = DirichletBC
    variable = disp_y
    boundary = 'left right top bottom'
    value = 0
  []
  [xfix]
    type = DirichletBC
    variable = disp_x
    boundary = 'left right top'
    value = 0
  []
[]

[Materials]
  [pressure]
    type = GenericFunctionMaterial
    prop_names = 'fracture_pressure'
    prop_values = 't'
  []
  [pfbulkmat]
    type = GenericConstantMaterial
    prop_names = 'gc_prop l visco'
    prop_values = '0.001 0.04 1e-6'
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
    function = 'gc_prop * l'
  []
  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1
    poissons_ratio = 0.2
  []
  [damage_stress]
    type = ComputeLinearElasticPFFractureStress
    c = c
    E_name = 'elastic_energy'
    D_name = 'degradation'
    F_name = 'local_fracture_energy'
    I_name = 'indicator_function'
    decomposition_type = none
    use_current_history_variable = true
    use_snes_vi_solver = true
  []
  [indicator_function]
    type = DerivativeParsedMaterial
    f_name = indicator_function
    args = 'c'
    function = 'c'
    derivative_order = 2
  []
  [degradation]
    type = DerivativeParsedMaterial
    f_name = degradation
    args = 'c'
    function = '(1.0-c)^2*(1.0 - eta) + eta'
    constant_names = 'eta'
    constant_expressions = '1e-8'
    derivative_order = 2
  []
  [local_fracture_energy]
    type = DerivativeParsedMaterial
    f_name = local_fracture_energy
    args = 'c'
    material_property_names = 'gc_prop l'
    function = 'c^2 * gc_prop / 2 / l'
    derivative_order = 2
  []
[]

[Postprocessors]
  [integrated_c]
    type = ElementIntegralVariablePostprocessor
    variable = 'c'
  []
[]

[Executioner]
  type = Transient

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -snes_type'
  petsc_options_value = 'lu       superlu_dist                  vinewtonrsls'

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8

  dt = 1e-3
  automatic_scaling = true

  picard_max_its = 100
  picard_rel_tol = 1e-6
  picard_abs_tol = 1e-8
  accept_on_max_picard_iteration = true
  end_time = ${end}
[]

[Outputs]
  print_linear_converged_reason = false
  print_nonlinear_converged_reason = false
  print_linear_residuals = false
  file_base = sneddon_a02
  exodus = true
  csv = true
[]
