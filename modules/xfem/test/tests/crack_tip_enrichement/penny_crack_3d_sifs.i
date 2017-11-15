[XFEM]
  geometric_cut_userobjects = 'circle_cut_uo'
  qrule = volfrac
  output_cut_plane = true
  use_crack_tip_enrichment = true
  crack_front_definition = crack_front
  enrichment_displacements = 'enrich1_x enrich2_x enrich3_x enrich4_x enrich1_y enrich2_y enrich3_y enrich4_y enrich1_z enrich2_z enrich3_z enrich4_z'
  displacements = 'disp_x disp_y disp_z'
  cut_off_boundary = all
  cut_off_radius = 0.08
[]

[UserObjects]
  [./circle_cut_uo]
    type = CircleCutUserObject
    cut_data = '0 0 0
                0.5 0 0
                0 0.5 0'
  [../]
  [./crack_front]
    type = CrackFrontDefinition
    crack_direction_method = CurvedCrackFront
    closed_loop = true
    crack_front_points = '   0.500000000000000                   0                   0
   0.492403876506104   0.086824088833465                   0
   0.469846310392954   0.171010071662834                   0
   0.433012701892219   0.250000000000000                   0
   0.383022221559489   0.321393804843270                   0
   0.321393804843270   0.383022221559489                   0
   0.250000000000000   0.433012701892219                   0
   0.171010071662834   0.469846310392954                   0
   0.086824088833465   0.492403876506104                   0
   0.000000000000000   0.500000000000000                   0
  -0.086824088833465   0.492403876506104                   0
  -0.171010071662834   0.469846310392954                   0
  -0.250000000000000   0.433012701892219                   0
  -0.321393804843270   0.383022221559489                   0
  -0.383022221559489   0.321393804843270                   0
  -0.433012701892219   0.250000000000000                   0
  -0.469846310392954   0.171010071662834                   0
  -0.492403876506104   0.086824088833465                   0
  -0.500000000000000   0.000000000000000                   0
  -0.492403876506104  -0.086824088833465                   0
  -0.469846310392954  -0.171010071662834                   0
  -0.433012701892219  -0.250000000000000                   0
  -0.383022221559489  -0.321393804843270                   0
  -0.321393804843270  -0.383022221559489                   0
  -0.250000000000000  -0.433012701892219                   0
  -0.171010071662835  -0.469846310392954                   0
  -0.086824088833465  -0.492403876506104                   0
  -0.000000000000000  -0.500000000000000                   0
   0.086824088833465  -0.492403876506104                   0
   0.171010071662834  -0.469846310392954                   0
   0.250000000000000  -0.433012701892219                   0
   0.321393804843270  -0.383022221559489                   0
   0.383022221559489  -0.321393804843270                   0
   0.433012701892219  -0.250000000000000                   0
   0.469846310392954  -0.171010071662835                   0
   0.492403876506104  -0.086824088833465                   0'
    [../]
[]

[DomainIntegral]
  integrals = 'Jintegral InteractionIntegralKI'
  convert_J_to_K = true
  displacements = 'disp_x disp_y disp_z'
  closed_loop = true
  crack_front_points = '   0.500000000000000                   0                   0
   0.492403876506104   0.086824088833465                   0
   0.469846310392954   0.171010071662834                   0
   0.433012701892219   0.250000000000000                   0
   0.383022221559489   0.321393804843270                   0
   0.321393804843270   0.383022221559489                   0
   0.250000000000000   0.433012701892219                   0
   0.171010071662834   0.469846310392954                   0
   0.086824088833465   0.492403876506104                   0
   0.000000000000000   0.500000000000000                   0
  -0.086824088833465   0.492403876506104                   0
  -0.171010071662834   0.469846310392954                   0
  -0.250000000000000   0.433012701892219                   0
  -0.321393804843270   0.383022221559489                   0
  -0.383022221559489   0.321393804843270                   0
  -0.433012701892219   0.250000000000000                   0
  -0.469846310392954   0.171010071662834                   0
  -0.492403876506104   0.086824088833465                   0
  -0.500000000000000   0.000000000000000                   0
  -0.492403876506104  -0.086824088833465                   0
  -0.469846310392954  -0.171010071662834                   0
  -0.433012701892219  -0.250000000000000                   0
  -0.383022221559489  -0.321393804843270                   0
  -0.321393804843270  -0.383022221559489                   0
  -0.250000000000000  -0.433012701892219                   0
  -0.171010071662835  -0.469846310392954                   0
  -0.086824088833465  -0.492403876506104                   0
  -0.000000000000000  -0.500000000000000                   0
   0.086824088833465  -0.492403876506104                   0
   0.171010071662834  -0.469846310392954                   0
   0.250000000000000  -0.433012701892219                   0
   0.321393804843270  -0.383022221559489                   0
   0.383022221559489  -0.321393804843270                   0
   0.433012701892219  -0.250000000000000                   0
   0.469846310392954  -0.171010071662835                   0
   0.492403876506104  -0.086824088833465                   0'
  crack_direction_method = CurvedCrackFront
  radius_inner = '0.075 0.1'
  radius_outer = '0.15 0.175'
  poissons_ratio = 0.3
  youngs_modulus = 1.0e6
  block = 0
  convert_J_to_K = true
[]


[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 51
  ny = 51
  nz = 51
  xmin = -1.0
  xmax = 1.0
  ymin = -1.0
  ymax = 1.0
  zmin = -1.0
  zmax = 1.0
  elem_type = HEX8
[]

[MeshModifiers]
  [./all_node]
    type = BoundingBoxNodeSet
    new_boundary = 'all'
    top_right = '1 1 1'
    bottom_left = '-1 -1 -1'
  [../]
[]

[Variables]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_z]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxVariables]
  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./SED]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./stress_xx]
  type = RankTwoAux
  rank_two_tensor = stress
  variable = stress_xx
  index_i = 0
  index_j = 0
  execute_on = timestep_end
  [../]
  [./stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
    execute_on = timestep_end
  [../]
  [./stress_zz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zz
    index_i = 2
    index_j = 2
    execute_on = timestep_end
  [../]
[]

[Kernels]
  [./TensorMechanics]
    displacements = 'disp_x disp_y disp_z'
    use_displaced_mesh = false
    volumetric_locking_correction = false
  [../]
[]

[BCs]
  [./top_z]
    type = Pressure
    variable = disp_z
    boundary = front
    component = 2
    factor = -1
  [../]
  [./bottom_x]
    type = PresetBC
    boundary = back
    variable = disp_x
    value = 0.0
  [../]
  [./bottom_y]
    type = PresetBC
    boundary = back
    variable = disp_y
    value = 0.0
  [../]
  [./bottom_z]
    type = PresetBC
    boundary = back
    variable = disp_z
    value = 0.0
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1e6
    poissons_ratio = 0.3
  [../]
  [./strain]
    type = ComputeCrackTipEnrichmentSmallStrain
    displacements = 'disp_x disp_y disp_z'
    crack_front_definition = crack_front
    enrichment_displacements = 'enrich1_x enrich2_x enrich3_x enrich4_x enrich1_y enrich2_y enrich3_y enrich4_y enrich1_z enrich2_z enrich3_z enrich4_z'
  [../]
  [./stress]
    type = ComputeLinearElasticStress
  [../]
  [./eshelby]
    type = EshelbyTensor
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

  solve_type = 'PJFNK'
#  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
#  petsc_options_value = 'lu     superlu_dist'

  petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type -pc_hypre_boomeramg_max_iter'
  petsc_options_value = '201                hypre    boomeramg      20'

  line_search = 'bt'

  [./Quadrature]
    type = GAUSS
    order = SECOND
  [../]

  # controls for linear iterations
  l_max_its = 25
  l_tol = 1e-3

  # controls for nonlinear iterations
  nl_max_its = 15
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-8

  # time control
  start_time = 0.0
  dt = 1.0
  end_time = 1.0
[]

[Outputs]
  exodus = true
  [./console]
    type = Console
    perf_log = true
    output_linear = true
  [../]
[]
