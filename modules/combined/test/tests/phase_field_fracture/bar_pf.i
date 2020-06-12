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
    bottom_left = '99.5 0 0'
    top_right = '100.5 1 0'
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
  [../]
  [./disp_y]
  [../]
[]

[Kernels]
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

[Materials]
  [./strain]
    type = ComputeSmallStrain
  [../]
  [./pfbulkmat]
    type = GenericConstantMaterial
    prop_names = 'gc_prop l visco'
    prop_values = '0.12 5 1.0e-10'
  [../]
  [./define_mobility]
    type = ParsedMaterial
    material_property_names = 'gc_prop visco'
    f_name = L
    function = '1.0'
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
    decomposition_type = none
    use_snes_vi_solver = true
    outputs = all
  [../]
  [./degradation]
    type = DerivativeParsedMaterial
    f_name = degradation
    args = 'c'
    function = '(1.0-c)^2/((1.0-c)^2+c*(1-0.5*c)*(4/3.14159/l*E*gc_prop/sigma^2))'
    material_property_names = 'gc_prop l'
    constant_names       = 'E sigma'
    constant_expressions = '30000 3'
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

  solve_type = PJFNK
  petsc_options_iname = '-pc_type  -snes_type'
  petsc_options_value = 'lu vinewtonrsls'

  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-12
  l_max_its = 10
  nl_max_its = 30

  dt = 1e-4
  end_time = 0.05
  dtmin = 1e-5
  automatic_scaling = true
[]

[Outputs]
  exodus = true
  csv = true
[]
