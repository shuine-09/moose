[Mesh]
  type = GeneratedMesh
  nx = 1
  ny = 1
  nz = 1
  dim = 3
  elem_type = HEX8
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
  [./ep_norm]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./e_zz]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./rho_i]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./rho_m]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./resid_z]
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
  [./TensorMechanics]
    displacements = 'disp_x disp_y disp_z'
    save_in_disp_z = resid_z
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
  [./rho_m]
    type = MaterialStdVectorAux
    variable = rho_m
    property = gss
    index = 14
    execute_on = timestep_end
    block = 0
  [../]
  [./rho_i]
    type = MaterialStdVectorAux
    variable = rho_i
    property = internal_var
    index = 14
    execute_on = timestep_end
    block = 0
  [../]
  [./ep_norm]
    type = RankTwoScalarAux
    variable = ep_norm
    rank_two_tensor = plastic_strain
    scalar_type = L2norm
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
  active = 'crysp'
  [./crysp]
    type = FiniteStrainCPDislocation
    block = 0
    disp_y = disp_y
    disp_x = disp_x
    gtol = 1e-5
    slip_sys_file_name = input_slip_sys.txt
    disp_z = disp_z
    C_ijkl = '279.06e3 119.6e3 119.6e3 279.06e3 119.6e3 279.06e3 79.731e3 79.731e3 79.731e3'
    nss = 48
    fill_method = symmetric9
    num_internal_var_nss = 4
    num_internal_var_scalar = 2
    intvar_read_type = slip_sys_file
    num_slip_sys_props = 4
    slip_incr_tol = 0.1
    rtol = 1e-6
    penalty_param = 1e-5
    burgers_length  = 0.248e-6 #mm
    disloc_line_dist = 1.5e-6  #mm
    jump_freq = 1e6
    active_enthal  = 1.2161e-18
    thermal_resist = 390.0
    temp = 323.0
    elastic_const_from_tensor = false
    maximum_substep_iteration = 9
  [../]
  [./elastic]
    type = FiniteStrainElasticMaterial
    block = 0
    disp_y = disp_y
    disp_x = disp_x
    disp_z = disp_z
    C_ijkl = '1.684e5 1.214e5 1.214e5 1.684e5 1.214e5 1.684e5 0.754e5 0.754e5 0.754e5'
    fill_method = symmetric9
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
  [./ep_norm]
    type = ElementAverageValue
    variable = ep_norm
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./resid_z]
    type = NodalSum
    variable = resid_z
    boundary = front
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
  solve_type = 'PJFNK'

  petsc_options_iname = -pc_hypre_type
  petsc_options_value = boomerang
  dtmax = 10.0
  nl_rel_tol = 1e-6
  nl_rel_step_tol = 1e-6
  end_time = 2.0
  dtmin = 1e-8
  num_steps = 5
  l_max_its = 60
  nl_max_its = 20

  [./TimeStepper]
   type = IterationAdaptiveDT
   reset_dt = true
   dt = 0.1
   time_t = '0.0 0.3'
   time_dt = '0.1 5e-3'
   growth_factor = 1.8
   cutback_factor = 0.5
  [../]
[]

[Outputs]
  file_base = crysp_substep_out
  exodus = true
  csv = true
  gnuplot = true
  [./console]
    type = Console
    perf_log = true
    output_linear = false
  [../]
[]

[Problem]
  use_legacy_uo_initialization = false
[]
