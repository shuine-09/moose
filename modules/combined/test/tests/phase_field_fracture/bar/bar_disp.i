#This input uses PhaseField-Nonconserved Action to add phase field fracture bulk rate kernels
[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 200
    ny = 1
    xmax = 200
    ymax = 1
  []
  [./center]
    type = BoundingBoxNodeSetGenerator
    input = 'gen'
    bottom_left = '-0.5 0 0'
    top_right = '0.5 1 0'
    new_boundary = 'center'
  [../]
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[MultiApps]
  [damage]
    type = TransientMultiApp
    input_files = 'bar_c.i'
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

[ICs]
  [./c]
    type = ConstantIC
    boundary = center
    value = 1e-3
    variable = c
  [../]
[]

# [Problem]
#   type = ReferenceResidualProblem
#       reference_vector = 'ref'
#         extra_tag_vectors = 'ref'
# []

[Modules]
  [./TensorMechanics]
    [./Master]
      [./mech]
        add_variables = true
        strain = SMALL
        additional_generate_output = 'stress_xx stress_zz'
        save_in = 'resid_x resid_y'
        #extra_vector_tags = 'ref'
      [../]
    [../]
  [../]
[]

[AuxVariables]
  [./c]
    order = FIRST
    family = LAGRANGE
  [../]
  [./resid_x]
  [../]
  [./resid_y]
  [../]
[]

[Kernels]
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
[]

[BCs]
  [./yfix]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  [../]
  [./xfix]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  [../]
  [./xdisp]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = right
    function = t
  [../]
[]

[Functions]
  [./fracture_pressure]
    type = ParsedFunction
    value = '0.4'
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
    prop_names = 'gc_prop visco'
    prop_values = '0.12 1.0'
  [../]
  [./pfbulkmat2]
    type = GenericConstantMaterial
    prop_names = 'l'
    prop_values = '20'
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
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 30000
    poissons_ratio = 0.2
  [../]
  [./elastic]
    type = ComputeLinearElasticPFFractureStress
    c = c
    E_name = 'elastic_energy'
    F_name = 'fracture_energy'
    D_name = 'degradation'
    decomposition_type = strain_vol_dev
    use_snes_vi_solver = true
    use_current_history_variable = true
    #outputs = all
  [../]
  [./degradation]
    type = DerivativeParsedMaterial
    f_name = degradation
    args = 'c'
    function = '(1.0-c)^2*(1-eta)/((1.0-c)^2+c*(1-0.5*c)*(4/3.14159/l*E*gc_prop/sigma^2))+eta'
    material_property_names = 'gc_prop l'
    constant_names       = 'E sigma eta'
    constant_expressions = '30000 3 1e-12'
    derivative_order = 2
  [../]
  [./fracture_energy]
    type = DerivativeParsedMaterial
    f_name = fracture_energy
    args = 'c'
    material_property_names = 'gc_prop l'
    function = 'gc_prop/l/3.14159*(2*c-c^2)'
    derivative_order = 2
  [../]
  [./fracture_driving_energy]
    type = DerivativeSumMaterial
    args = c
    sum_materials = 'elastic_energy fracture_energy'
    derivative_order = 2
    f_name = F
    #outputs = all
  [../]
[]

[Postprocessors]
  [./resid_x_left]
    type = NodalSum
    variable = resid_x
    boundary = left
  [../]

  [./resid_x]
    type = NodalSum
    variable = resid_x
    boundary = right
  [../]
  [./resid_y]
    type = NodalSum
    variable = resid_y
    boundary = right
  [../]
  [./disp_x]
    type = NodalMaxValue
    variable = disp_x
    boundary = right
  [../]
  [./pressure]
    type = ElementIntegralVariablePostprocessor2
    displacements = 'disp_x disp_y'
    variable = c
  [../]
  [./max_d]
    type = ElementExtremeValue
    variable = c
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount '
  petsc_options_value = 'lu  NONZERO 1e-10'

  line_search = none

  nl_rel_tol = 1e-7
  nl_abs_tol = 1e-9
  l_max_its = 10
  nl_max_its = 10

  picard_max_its = 1000
  picard_rel_tol = 1e-7
  picard_abs_tol = 1e-7
  accept_on_max_picard_iteration = true

  #abort_on_solve_fail = true

  dt = 1e-4
  end_time = 0.08
  dtmin = 1e-14
  automatic_scaling = true
[]

[Outputs]
  file_base = bar_p00_l50_3
#exodus = true
  csv = true
[]
