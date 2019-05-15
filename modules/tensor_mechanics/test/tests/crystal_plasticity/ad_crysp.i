[Mesh]
  type = GeneratedMesh
  dim = 3
  elem_type = HEX8
  nx = 4
  ny = 4
  nz = 4
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Variables]
  [./disp_x]
    block = 0
  [../]
  [./disp_y]
    block = 0
  [../]
  [./disp_z]
    block = 0
  [../]
[]

[AuxVariables]
  [./stress_zz]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./fp_zz]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./rotout]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./e_zz]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
[]

[Functions]
  [./tdisp]
    type = ParsedFunction
    value = 0.01*t
  [../]
[]

[Kernels]
  [./stress_x]
    type = ADStressDivergenceTensors
    component = 0
    variable = disp_x
    use_displaced_mesh = true
  [../]
  [./stress_y]
    type = ADStressDivergenceTensors
    component = 1
    variable = disp_y
    use_displaced_mesh = true
  [../]
  [./stress_z]
    type = ADStressDivergenceTensors
    component = 2
    variable = disp_z
    use_displaced_mesh = true
  [../]
[]

[AuxKernels]
  [./stress_zz]
    type = RankTwoAux
    variable = stress_zz
    rank_two_tensor = stress
    index_j = 2
    index_i = 2
    execute_on = timestep_end
    block = 0
  [../]
  [./fp_zz]
    type = RankTwoAux
    variable = fp_zz
    rank_two_tensor = fp
    index_j = 2
    index_i = 2
    execute_on = timestep_end
    block = 0
  [../]
  [./e_zz]
    type = RankTwoAux
    variable = e_zz
    rank_two_tensor = lage
    index_j = 2
    index_i = 2
    execute_on = timestep_end
    block = 0
  [../]
[]

[BCs]
  [./symmy]
    type = PresetBC
    variable = disp_y
    boundary = bottom
    value = 0
  [../]
  [./symmx]
    type = PresetBC
    variable = disp_x
    boundary = left
    value = 0
  [../]
  [./symmz]
    type = PresetBC
    variable = disp_z
    boundary = back
    value = 0
  [../]
  [./tdisp]
    type = FunctionPresetBC
    variable = disp_z
    boundary = front
    function = tdisp
  [../]
[]

[Materials]
  [./crysp]
    type = ADFiniteStrainCrystalPlasticity
    block = 0
    gtol = 1e-2
    slip_sys_file_name = input_slip_sys.txt
    nss = 12
    num_slip_sys_flowrate_props = 2 #Number of properties in a slip system
    flowprops = '1 4 0.001 0.1 5 8 0.001 0.1 9 12 0.001 0.1'
    hprops = '1.0 541.5 60.8 109.8 2.5'
    gprops = '1 4 60.8 5 8 60.8 9 12 60.8'
    maximum_substep_iteration = 1
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensorCP
    block = 0
    C_ijkl = '1.684e5 1.214e5 1.214e5 1.684e5 1.214e5 1.684e5 0.754e5 0.754e5 0.754e5'
    fill_method = symmetric9
  [../]
  [./strain]
    type = ADComputeGreenLagrangeStrain
    block = 0
  [../]
[]

[Postprocessors]
  [./stress_zz]
    type = ElementAverageValue
    variable = stress_zz
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./fp_zz]
    type = ElementAverageValue
    variable = fp_zz
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./e_zz]
    type = ElementAverageValue
    variable = e_zz
    block = 'ANY_BLOCK_ID 0'
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

  #Preconditioned JFNK (default)
  solve_type = 'NEWTON'

  petsc_options_iname = '-pc_type -ksp_type'
  petsc_options_value = 'lu preonly'

  line_search = 'none'

  nl_rel_tol = 1e-11
  nl_abs_tol = 1e-11

  nl_max_its = 10

  dt = 0.1
  dtmax = 10.0
  dtmin = 0.01

  num_steps = 10
[]

[Outputs]
  file_base = out
  exodus = true
  perf_graph = true
[]
