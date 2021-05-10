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

[Variables]
  [./c]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[ICs]
  [./c]
    type = ConstantIC
    boundary = center
    value = 1e-3
    variable = c
  [../]
[]

[AuxVariables]
  [./bounds_dummy]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
#  [./pf_c]
#    type = PhaseFieldPressurizedFractureDamage
#    variable = c
#    displacements = 'disp_x disp_y'
#  [../]
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
  [strain]
    type = ComputeSmallStrain
  []
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
    decomposition_type = none
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
    outputs = all
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
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

[Executioner]
  type = Transient

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -snes_type'
  petsc_options_value = 'lu       superlu_dist                  vinewtonrsls'

  line_search = none

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  l_max_its = 10
  nl_max_its = 100

  dt = 1e-3
  end_time = 0.08
  dtmin = 1e-10
  automatic_scaling = true
[]

[Outputs]
  print_linear_converged_reason = false
  print_nonlinear_converged_reason = false
  print_linear_residuals = true
[]
