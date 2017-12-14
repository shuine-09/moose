[Mesh]
  type = FileMesh
  file = sent7.e
  uniform_refine = 0
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Modules]
  [./TensorMechanics]
    [./Master]
      [./mech]
        add_variables = true
        strain = SMALL
        additional_generate_output = stress_yy
        decomposition_method = EigenSolution
        save_in = 'resid_x resid_y'
      [../]
    [../]
  [../]
[]

[Variables]
  [./c]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxVariables]
  [./resid_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./resid_y]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./ac_pf]
    type = AllenCahnPFFracture
    l_name = l
    visco_name = visco
    gc = gc_prop
    displacements = 'disp_x disp_y'
    F_name = E_el
    variable = c
  [../]
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

[Materials]
  [./pfbulkmat]
    type = GenericConstantMaterial
    prop_names = 'gc_prop l visco'
    prop_values = '1e-3 0.01 1e-3'
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
    kdamage = 1e-03
  [../]
  #[./elasticity_tensor]
  #  type = ComputeElasticityTensor
  #  C_ijkl = '120.0 80.0'
  #  fill_method = symmetric_isotropic
  #[../]
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    C_ijkl = '1.684e2 1.214e2 1.214e2 1.684e2 1.214e2 1.684e2 0.754e2 0.754e2 0.754e2'
    fill_method = symmetric9
    euler_angle_1 = 0
    euler_angle_2 = 0
    euler_angle_3 = 0
  [../]
[]

[BCs]
  [./xdisp]
    type = FunctionPresetBC
    variable = disp_y
    boundary = 2
    function = 't'
  [../]
  [./xfix]
    type = PresetBC
    variable = disp_y
    boundary = 1
    value = 0
  [../]
  [./yfix]
    type = PresetBC
    variable = disp_x
    boundary = '1 2'
    value = 0
  [../]
[]

[Preconditioning]
  active = 'smp'
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Dampers]
  [./bounding_value_damp]
    type = BoundingValueElementDamper
    min_value = -0.01
    max_value = 1.02
    variable = c
  [../]
[]

[Postprocessors]
  [./disp_x_top]
    type = SideAverageValue
    variable = disp_y
    boundary = 2
  [../]
  [./reaction_force_x]
    type = NodalSum
    variable = resid_y
    boundary = 2
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

  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-6
  l_tol = 1e-4
  l_max_its = 50
  nl_max_its = 50

  dt = 2.5e-5
  num_steps = 1000
[]

[Outputs]
#  file_base = fd_test
  file_base = test_r0_mode1
  exodus = true
  csv = true
  print_perf_log = true
[]
