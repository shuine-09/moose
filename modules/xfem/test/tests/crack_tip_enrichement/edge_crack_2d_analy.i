[XFEM]
  geometric_cut_userobjects = 'line_seg_cut_uo'
  qrule = volfrac
  output_cut_plane = true
  use_crack_tip_enrichment = true
  crack_front_definition = crack_tip
  enrichment_displacements = 'enrich1_x enrich2_x enrich3_x enrich4_x enrich1_y enrich2_y enrich3_y enrich4_y'
  cut_off_boundary = all
  cut_off_radius = 0.1
[]

[UserObjects]
  [./line_seg_cut_uo]
    type = LineSegmentCutUserObject
    cut_data = '0.0 1.0 0.5 1.0'
    time_start_cut = 0.0
    time_end_cut = 0.0
  [../]
  [./crack_tip]
    type = CrackFrontDefinition
    crack_direction_method = CrackDirectionVector
    crack_front_points = '0.5 1.0 0'
    crack_direction_vector = '1 0 0'
    2d = true
    axis_2d = 2
  [../]
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 41
  ny = 81
  xmin = 0.0
  xmax = 1.0
  ymin = 0.0
  ymax = 2.0
  elem_type = QUAD4
  displacements = 'disp_x disp_y'
[]

[AuxVariables]
  [./SED]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./SED]
    type = MaterialRealAux
    variable = SED
    property = strain_energy_density
    execute_on = timestep_end
  [../]
[]

[DomainIntegral]
  integrals = JIntegral
  convert_J_to_K = true
  crack_front_points = '0.5 1.0 0'
  crack_direction_method = CrackDirectionVector
  crack_direction_vector = '1 0 0'
  2d = true
  axis_2d = 2
  radius_inner = '0.05 0.06 0.07 0.08 0.09 0.1 0.11 0.12'
  radius_outer = '0.1 0.12 0.14 0.16 0.18 0.2 0.22 0.24'
  #radius_inner = '0.05'
  #radius_outer = '0.1'
  youngs_modulus = 1e6
  poissons_ratio = 0.3
[]

[MeshModifiers]
  [./all_node]
    type = BoundingBoxNodeSet
    new_boundary = 'all'
    top_right = '1 2 0'
    bottom_left = '0 0 0'
  [../]
  [./right_bottom_node]
    type = AddExtraNodeset
    new_boundary = 'right_bottom_node'
    coord = '1.0 0.0'
  [../]
  [./right_top_node]
    type = AddExtraNodeset
    new_boundary = 'right_top_node'
    coord = '1.0 2.0'
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
[]

[AuxVariables]
 [./saved_x]
  [../]
  [./saved_y]
  [../]
  [./stress_xx]      # stress aux variables are defined for output; this is a way to get integration point variables to the output file
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vonmises]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  [./TensorMechanics]
    displacements = 'disp_x disp_y'
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
  [./stress_xy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xy
    index_i = 0
    index_j = 1
    execute_on = timestep_end
  [../]
  [./vonmises]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = vonmises
    scalar_type = vonmisesStress
    execute_on = timestep_end
  [../]
[]

[BCs]
  [./top_y]
    type = Pressure
    variable = disp_y
    boundary = top
    component = 1
    factor = -1
  [../]
  [./bottom_y]
    type = Pressure
    variable = disp_y
    boundary = bottom
    component = 1
    factor = -1
  [../]
  [./fix_y]
    type = PresetBC
    boundary = right_bottom_node
    variable = disp_y
    value = 0.0
  [../]
  [./fix_x]
    type = PresetBC
    boundary = right_bottom_node
    variable = disp_x
    value =  0.0
  [../]
  [./fix_x2]
    type = PresetBC
    boundary = right_top_node
    variable = disp_x
    value =  0.0
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
    displacements = 'disp_x disp_y'
    crack_front_definition = crack_tip
    enrichment_displacements = 'enrich1_x enrich2_x enrich3_x enrich4_x enrich1_y enrich2_y enrich3_y enrich4_y'
  [../]
  [./stress]
    type = ComputeLinearElasticStress
  [../]
  [./eshelby]
    type = EshelbyTensor
    displacements = 'disp_x disp_y'
  [../]
[]

[Postprocessors]
  [./error]
    type = MaterialTensorIntegralXFEM
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
    execute_on = timestep_end
 [../]
 [./error2]
   type = MaterialTensorIntegralXFEM
   rank_two_tensor = stress
   index_i = 1
   index_j = 1
   execute_on = timestep_end
   normalized_error = true
[../]
[]

[Executioner]
  type = Transient

  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu     superlu_dist'

  line_search = 'none'

  [./Quadrature]
    type = GAUSS
    order = SIXTH
  [../]

  [./Predictor]
    type = SimplePredictor
    scale = 1.0
  [../]

# controls for linear iterations
  l_max_its = 100
  l_tol = 1e-4

# controls for nonlinear iterations
  nl_max_its = 100
  nl_rel_tol = 1e-9 #11
  nl_abs_tol = 1e-8 #12

# time control
  start_time = 0.0
  dt = 1.0
  end_time = 1.0
  dtmin = 1.0
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]


[Outputs]
  exodus = true
  [./console]
    type = Console
    perf_log = true
    output_linear = true
  [../]
[]
