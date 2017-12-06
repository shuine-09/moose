[Mesh]
  type = FileMesh
  file = sent6.e
  uniform_refine = 0
[]
#[Mesh]
#  type = GeneratedMesh
#  nx = 200
#  ny = 80
#  xmax = 1.0
#  xmin = 0.0
#  ymax = 0.4
#  ymin = 0.0
#  dim = 2
#[]

[GlobalParams]
  displacements = 'disp_x disp_y'
  block = 1
[]

[AuxVariables]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Modules]
  [./PhaseField]
    [./Nonconserved]
      [./sub_c]
        free_energy = E_el
        kappa = kappa_op
        mobility = L
      [../]
    [../]
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
    c = sub_c
    F_name = E_el
    kdamage = 1.0e-4
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    C_ijkl = '120.0 80.0'
    fill_method = symmetric_isotropic
  [../]
  [./strain]
    type = ComputeSmallStrain
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
  #petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
  #petsc_options_value = 'asm      31                  preonly       lu           1'

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu     superlu_dist'

  line_search = default

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  l_tol = 1e-4
  l_max_its = 50
  nl_max_its = 30

  dt = 1.0e-4
  num_steps = 10000
[]

#[Dampers]
#  [./bounding_value_damp]
#    type = BoundingValueElementDamper
#    min_value = 0
#    max_value = 1.04
#    variable = sub_c
#  [../]
#[]

[Outputs]
  print_linear_residuals = true
  csv = true
  exodus = true
[]
