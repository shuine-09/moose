[Mesh]
  type = FileMesh
  file = L_3D.e
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Modules]
  [./TensorMechanics]
    [./Master]
      [./mech]
        add_variables = true
        strain = SMALL
        additional_generate_output = stress_yy
        decomposition_method = EigenSolution
        save_in = 'resid_x resid_y resid_z'
      [../]
    [../]
  [../]
[]

#[XFEM]
#  geometric_cut_userobjects = 'line_seg_cut_uo'
#  qrule = volfrac
#  output_cut_plane = true
#[]
#
#[UserObjects]
#  [./line_seg_cut_uo]
#    type = LineSegmentCutUserObject
#    cut_data = '-4.99 -4.0 -4.99 -2.5'
#    time_start_cut = 0.0
#    time_end_cut = 0.0
#  [../]
#[]

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
  [./resid_z]
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
    displacements = 'disp_x disp_y disp_z'
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
  [./solid_z]
    type = PhaseFieldFractureMechanicsOffDiag
    variable = disp_z
    component = 2
    c = c
  [../]
[]

[Materials]
  [./pfbulkmat]
    type = GenericConstantMaterial
    prop_names = 'gc_prop l visco'
    prop_values = '9.5e-2 12.0 1e-3'
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
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    C_ijkl = '2.58e4 1.09e4'
    fill_method = symmetric_isotropic
  [../]
[]

[BCs]
  [./ydisp]
    type = FunctionPresetBC
    variable = disp_y
    boundary = load
    function = disp_y_t
  [../]
  [./bottom_x]
    type = PresetBC
    variable = disp_x
    boundary = bottom
    value = 0
  [../]
  [./bottom_y]
    type = PresetBC
    variable = disp_y
    boundary = bottom
    value = 0
  [../]
  [./bottom_z]
    type = PresetBC
    variable = disp_z
    boundary = bottom
    value = 0
  [../]
[]

[Functions]
  [./disp_y_t]
    type = ParsedFunction
    value = t
  [../]
[]

[Postprocessors]
  [./reaction_x]
    type = NodalSum
    variable = resid_x
    boundary = load
  [../]
  [./reaction_y]
    type = NodalSum
    variable = resid_y
    boundary = load
  [../]
  [./reaction_z]
    type = NodalSum
    variable = resid_z
    boundary = load
  [../]
[]

[Preconditioning]
  active = 'smp'
  [./smp]
    type = SMP
    full = true
  [../]
[]

#[Dampers]
#  [./bounding_value_damp]
#    type = BoundingValueElementDamper
#    min_value = -0.01
#    max_value = 1.02
#    variable = c
#  [../]
#[]


[Executioner]
  type = Transient

  solve_type = PJFNK
  #petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
  #petsc_options_value = 'asm      31                  preonly       lu           1'

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu       superlu_dist'
  nl_abs_tol = 1e-5
  nl_rel_tol = 1e-5
  l_max_its = 100
  nl_max_its = 30

  dt = 1e-3
  dtmin = 1e-10
  start_time = 0.0
  end_time = 1.0
[]

[Outputs]
  file_base = l_shape
  exodus = true
  csv = true
  print_perf_log = true
[]
