[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 200
    ny = 200
    xmin = 0
    xmax = 4
    ymin = 0
    ymax = 4
  []
[]

[ICs]
  # [./ini_c]
  #   type = PFFractureLineIC
  #   variable = c
  #   l = 0.02
  #   point_1 = '1.8 2 0'
  #   point_2 = '2.2 2 0'
  # [../]
  [./ini_c]
    type = BrittleDamageIC
    variable = c
    l = 0.08
    d0 = 1.0
    x1 = '1.8 1.8'
    x2 = '2.2 2.2'
    y1 = '1.5 2.5'
    y2 = '1.5 2.5'
    z1 = '0 0'
    z2 = '0 0'
  [../]
[]

# [Mesh]
#   type = FileMesh
#   file = visualize.e
#   parallel_type = replicated
# []
#
# # [ICs]
# #   [./c]
# #     type = FunctionIC
# #     function = initial_c
# #     variable = c
# #   [../]
# # []
# #
# [Functions]
#   [./initial_c]
#     type = SolutionFunction
#     solution = ex_soln
#   [../]
# []
#
# [UserObjects]
#   [./ex_soln]
#     type = SolutionUserObject
#     system_variables = d
#     mesh = visualize.e
#     timestep = 2
#   [../]
# []

[Modules]
  [./TensorMechanics]
    [./Master]
      [./mech]
        add_variables = true
        strain = SMALL
        additional_generate_output = 'stress_yy'
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
    order = FIRST
    family = MONOMIAL
  [../]
  [./unmodified_stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./elastic_energy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  # [./c]
  # [../]
[]

[AuxKernels]
  [./cod]
    type = PFCrackOpeningDisplacement
    displacements = 'disp_x disp_y'
    c = c
    variable = cod
  [../]
  [./unmodified_stress_yy]
    type = MaterialRankTwoTensorAux
    i = 1
    j = 1
    property = unmodified_stress
    variable = unmodified_stress_yy
  [../]
  [./elastic_energy]
    type = MaterialRealAux
    property = elastic_energy
    variable = elastic_energy
  [../]
  # [./c]
  #   type = FunctionAux
  #   function = initial_c
  #   variable = c
  # [../]
[]

[Kernels]
  # [./solid_x]
  #   type = PhaseFieldFractureMechanicsOffDiag
  #   variable = disp_x
  #   component = 0
  #   c = c
  # [../]
  # [./solid_y]
  #   type = PhaseFieldFractureMechanicsOffDiag
  #   variable = disp_y
  #   component = 1
  #   c = c
  # [../]
#   [./diff]
#   type = Diffusion
#   displacements = 'disp_x disp_y'
#   variable = c
# [../]
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
  [./v_bounds]
    type = PFFractureBoundsAux
    variable = bounds_dummy
    bounded_variable = c
    upper = 1.0
    lower = 0
    lower_var = c
  [../]
[]

[Functions]
  [./fracture_pressure]
    type = ParsedFunction
    value = '(0.1+t*0.1)'
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
    prop_values = '1 0.08 1e-4'
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
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    C_ijkl = '1 0.2'
    fill_method = symmetric_isotropic_E_nu
  [../]
  [./damage_stress]
    type = ComputeLinearElasticPFFractureStress
    c = c
    E_name = 'elastic_energy'
    D_name = 'degradation'
    I_name = 'indicator_function'
    F_name = 'local_fracture_energy'
    decomposition_type = strain_spectral
    use_current_history_variable = false
    #decomposition_type = none
    use_snes_vi_solver = true
  [../]
  [./degradation]
    type = DerivativeParsedMaterial
    f_name = degradation
    args = 'c'
    function = '(1.0-c)^2*(1.0 - eta) + eta'
    constant_names       = 'eta'
    constant_expressions = '1.0e-6'
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

[Executioner]
  type = Transient

  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu       superlu_dist'

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-6
  l_max_its = 10
  nl_max_its = 25
  line_search = 'none'

  dt = 0.5
  #num_steps = 2
  end_time = 30
  automatic_scaling = true
  # [./Adaptivity]
  #  initial_adaptivity = 2
  #   refine_fraction = 0.99
  #   coarsen_fraction = 0.0
  #   max_h_level = 2
  # [../]
[]

[Outputs]
  exodus = true
  csv = true
[]
