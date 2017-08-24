[GlobalParams]
  order = FIRST
  family = LAGRANGE
  crack_front_definition = crack_front
  enrichment_displacement_var = 'enrich1_x enrich1_y enrich1_z enrich2_x enrich2_y enrich2_z enrich3_x enrich3_y enrich3_z enrich4_x enrich4_y enrich4_z'
  enrichment_displacement = 'enrich1_x enrich1_y enrich1_z enrich2_x enrich2_y enrich2_z enrich3_x enrich3_y enrich3_z enrich4_x enrich4_y enrich4_z'
  cut_off_radius = 0.15
[]

[XFEM]
  geometric_cut_userobjects = 'circle_cut_uo'
  qrule = volfrac
  output_cut_plane = true
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
    crack_front_points = '0.500000000000000                   0                   0
                          0.475528258147577   0.154508497187474                   0
                          0.404508497187474   0.293892626146237                   0
                          0.293892626146237   0.404508497187474                   0
                          0.154508497187474   0.475528258147577                   0
                          0.000000000000000   0.500000000000000                   0
                         -0.154508497187474   0.475528258147577                   0
                         -0.293892626146237   0.404508497187474                   0
                         -0.404508497187474   0.293892626146237                   0
                         -0.475528258147577   0.154508497187474                   0
                         -0.500000000000000   0.000000000000000                   0
                         -0.475528258147577  -0.154508497187473                   0
                         -0.404508497187474  -0.293892626146237                   0
                         -0.293892626146237  -0.404508497187474                   0
                         -0.154508497187474  -0.475528258147577                   0
                         -0.000000000000000  -0.500000000000000                   0
                          0.154508497187474  -0.475528258147577                   0
                          0.293892626146236  -0.404508497187474                   0
                          0.404508497187474  -0.293892626146237                   0
                          0.475528258147577  -0.154508497187474                   0'
  [../]
[]

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 11
  ny = 11
  nz = 7
  xmin = -1.0
  xmax = 1.0
  ymin = -1.0
  ymax = 1.0
  zmin = -0.5
  zmax = 0.5
  elem_type = HEX8
  displacements = 'disp_x disp_y disp_z'
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
  [../]
  [./disp_y]
  [../]
  [./disp_z]
  [../]
  [./enrich1_x]
  [../]
  [./enrich1_y]
  [../]
  [./enrich1_z]
  [../]
  [./enrich2_x]
  [../]
  [./enrich2_y]
  [../]
  [./enrich2_z]
  [../]
  [./enrich3_x]
  [../]
  [./enrich3_y]
  [../]
  [./enrich3_z]
  [../]
  [./enrich4_x]
  [../]
  [./enrich4_y]
  [../]
  [./enrich4_z]
  [../]
  []

[AuxVariables]
  [./stress_xx]      # stress aux variables are defined for output; this is a way to get integration point variables to the output file
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
[]

[AuxKernels]
  [./stress_xx]               # computes stress components for output
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
  [./enrich1_x]
    type = CrackTipEnrichmentStressDivergenceTensors
    variable = enrich1_x
    component = 0
    enrichment_component = 0
  [../]
  [./enrich1_y]
    type = CrackTipEnrichmentStressDivergenceTensors
    variable = enrich1_y
    component = 1
    enrichment_component = 0
  [../]
  [./enrich1_z]
    type = CrackTipEnrichmentStressDivergenceTensors
    variable = enrich1_z
    component = 2
    enrichment_component = 0
  [../]
  [./enrich2_x]
    type = CrackTipEnrichmentStressDivergenceTensors
    variable = enrich2_x
    component = 0
    enrichment_component = 1
  [../]
  [./enrich2_y]
    type = CrackTipEnrichmentStressDivergenceTensors
    variable = enrich2_y
    component = 1
    enrichment_component = 1
  [../]
  [./enrich2_z]
    type = CrackTipEnrichmentStressDivergenceTensors
    variable = enrich2_z
    component = 2
    enrichment_component = 1
  [../]
  [./enrich3_x]
    type = CrackTipEnrichmentStressDivergenceTensors
    variable = enrich3_x
    component = 0
    enrichment_component = 2
  [../]
  [./enrich3_y]
    type = CrackTipEnrichmentStressDivergenceTensors
    variable = enrich3_y
    component = 1
    enrichment_component = 2
  [../]
  [./enrich3_z]
    type = CrackTipEnrichmentStressDivergenceTensors
    variable = enrich3_z
    component = 2
    enrichment_component = 2
  [../]
  [./enrich4_x]
    type = CrackTipEnrichmentStressDivergenceTensors
    variable = enrich4_x
    component = 0
    enrichment_component = 3
  [../]
  [./enrich4_y]
    type = CrackTipEnrichmentStressDivergenceTensors
    variable = enrich4_y
    component = 1
    enrichment_component = 3
  [../]
  [./enrich4_z]
    type = CrackTipEnrichmentStressDivergenceTensors
    variable = enrich4_z
    component = 2
    enrichment_component = 3
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
  [./enrich1x_bc]
    type = CrackTipEnrichmentCutOffBC
    boundary = all
    variable = enrich1_x
    value =  0
  [../]
  [./enrich1y_bc]
    type = CrackTipEnrichmentCutOffBC
    boundary = all
    variable = enrich1_y
    value =  0
  [../]
  [./enrich1z_bc]
    type = CrackTipEnrichmentCutOffBC
    boundary = all
    variable = enrich1_z
    value =  0
  [../]
  [./enrich2x_bc]
    type = CrackTipEnrichmentCutOffBC
    boundary = all
    variable = enrich2_x
    value =  0
  [../]
  [./enrich2y_bc]
    type = CrackTipEnrichmentCutOffBC
    boundary = all
    variable = enrich2_y
    value =  0
  [../]
  [./enrich2z_bc]
     type = CrackTipEnrichmentCutOffBC
     boundary = all
     variable = enrich2_z
     value =  0
   [../]
  [./enrich3x_bc]
    type = CrackTipEnrichmentCutOffBC
    boundary = all
    variable = enrich3_x
    value =  0
  [../]
  [./enrich3y_bc]
    type = CrackTipEnrichmentCutOffBC
    boundary = all
    variable = enrich3_y
    value =  0
  [../]
  [./enrich3z_bc]
    type = CrackTipEnrichmentCutOffBC
    boundary = all
    variable = enrich3_z
    value =  0
  [../]
  [./enrich4x_bc]
    type = CrackTipEnrichmentCutOffBC
    boundary = all
    variable = enrich4_x
    value =  0
  [../]
  [./enrich4y_bc]
    type = CrackTipEnrichmentCutOffBC
    boundary = all
    variable = enrich4_y
    value =  0
  [../]
  [./enrich4z_bc]
    type = CrackTipEnrichmentCutOffBC
    boundary = all
    variable = enrich4_z
    value =  0
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
  [../]
  [./stress]
    type = ComputeLinearElasticStress
  [../]
[]

[Executioner]
  type = Transient

  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu     superlu_dist'

  line_search = 'none'

  #[./Predictor]
  #  type = SimplePredictor
  #  scale = 1.0
  #[../]

  [./Quadrature]
    type = GAUSS
    order = SIXTH
  [../]

# controls for linear iterations
  l_max_its = 20
  l_tol = 1e-4

# controls for nonlinear iterations
  nl_max_its = 15
  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-6

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
