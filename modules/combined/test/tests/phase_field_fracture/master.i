[Mesh]
  type = FileMesh
  file = sent6.e
  uniform_refine = 0
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
  block = 1
[]

[AuxVariables]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./c]
    order = FIRST
    family = LAGRANGE
  [../]
[]


[Modules]
  [./TensorMechanics]
    [./Master]
      [./mech]
        add_variables = true
        strain = SMALL
        incremental = false
        decomposition_method = EigenSolution
      [../]
    [../]
  [../]
[]

[Kernels]
  [./solid_x]
    type = PhaseFieldFractureMechanicsOffDiag
    variable = disp_x
    component = 0
    c = c
  [../]
  [./solid_y]
    type = PhaseFieldFractureMechanicsOffDiag
    variable = disp_y
    component = 1
    c = c
  [../]
[]

#[Functions]
#  [./pull_up_and_down]
#    type = ParsedFunction
#    value = 'if(t < 0.002, 5 * t, 1 * (t-0.002) + 0.01)'
#  [../]
#[]

[AuxKernels]
  [./stress_yy]
    type = RankTwoAux
    variable = stress_yy
    rank_two_tensor = stress
    index_j = 1
    index_i = 1
    execute_on = timestep_end
  [../]
[]

[BCs]
  [./xdisp]
    type = FunctionPresetBC
    variable = disp_x
    boundary = 2
    function = 't'
    #function = pull_up_and_down
  [../]
  [./xfix]
    type = PresetBC
    variable = disp_x
    boundary = 1
    value = 0
  [../]
  [./yfix]
    type = PresetBC
    variable = disp_y
    boundary = '1 2'
    value = 0
  [../]
[]

[Materials]
  [./pfbulkmat]
    type = GenericConstantMaterial
    prop_names = 'gc_prop l visco'
    prop_values = '1e-3 0.0075 1e-3'
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
  [./elastic]
    type = ComputeLinearElasticPFFractureStress
    c = c
    F_name = E_el
    kdamage = 1.0e-4
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    C_ijkl = '120.0 80.0'
    fill_method = symmetric_isotropic
  [../]
[]

[Postprocessors]
  [./disp_y_top]
    type = SideAverageValue
    variable = disp_y
    boundary = 2
  [../]
  [./stressyy]
    type = ElementAverageValue
    variable = stress_yy
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[MultiApps]
  [./sub]
    type = TransientMultiApp
    app_type = CombinedApp
    execute_on = timestep_end
    positions = '0 0 0'
    input_files = sub.i
  [../]
[]

[Transfers]
  [./from_sub]
    type = MultiAppMeshFunctionTransfer
    direction = from_multiapp
    multi_app = sub
    source_variable = sub_c
    variable = c
  [../]
  [./to_sub_x]
    type = MultiAppMeshFunctionTransfer
    direction = to_multiapp
    multi_app = sub
    source_variable = disp_x
    variable = disp_x
  [../]
  [./to_sub_y]
    type = MultiAppMeshFunctionTransfer
    direction = to_multiapp
    multi_app = sub
    source_variable = disp_y
    variable = disp_y
  [../]
[]

[Executioner]
  type = Transient

  solve_type = PJFNK
  #petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
  #petsc_options_value = 'asm      31                  preonly       lu           1'

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu     superlu_dist'

  line_search = none

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6
  l_tol = 1e-4
  l_max_its = 50
  nl_max_its = 30

  dt = 1.0e-5
  num_steps = 10000
[]

[Outputs]
  print_linear_residuals = true
  csv = true
  exodus = true
[]
