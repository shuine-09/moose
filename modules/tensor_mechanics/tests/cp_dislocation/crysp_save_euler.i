[Mesh]
  type = GeneratedMesh
  dim = 2
  xmin = 0.0
  xmax = 1.0
  ymin = 0.0
  ymax = 1.0
  nx = 2
  ny = 2
  elem_type = QUAD4

  displacements = 'disp_x disp_y'
[]

[Variables]
  [./disp_x]
    block = 0
  [../]
  [./disp_y]
    block = 0
  [../]
[]

[AuxVariables]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./fp_yy]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./e_yy]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./rho_m]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./rho_i]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./euler1]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./euler2]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./euler3]
    order = CONSTANT
    family = MONOMIAL
    block = 0
  [../]
  [./resid_y]
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
    displacements = 'disp_x disp_y'
    save_in_disp_y = resid_y
  [../]
[]

[AuxKernels]
  [./stress_yy]
    type = RankTwoAux
    variable = stress_yy
    rank_two_tensor = stress
    index_j = 1
    index_i = 1
    execute_on = timestep_end
    block = 0
  [../]
  [./fp_yy]
    type = RankTwoAux
    variable = fp_yy
    rank_two_tensor = fp
    index_j = 1
    index_i = 1
    execute_on = timestep_end
    block = 0
  [../]
  [./e_yy]
    type = RankTwoAux
    variable = e_yy
    rank_two_tensor = lage
    index_j = 1
    index_i = 1
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
  [./euler1]
    type = MaterialStdVectorAux
    variable = euler1
    property = euler_ang
    index = 0
    execute_on = timestep_end
    block = 0
  [../]
  [./euler2]
    type = MaterialStdVectorAux
    variable = euler2
    property = euler_ang
    index = 1
    execute_on = timestep_end
    block = 0
  [../]
  [./euler3]
    type = MaterialStdVectorAux
    variable = euler3
    property = euler_ang
    index = 2
    execute_on = timestep_end
    block = 0
  [../]
[]

[BCs]
  [./yeq0]
    type = PresetBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]
  [./xeq0]
    type = PresetBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]
  [./tdisp]
    type = FunctionPresetBC
    variable = disp_y
    boundary = top
    function = tdisp
  [../]
[]

[Materials]
  [./crysp]
    type = FiniteStrainCPDislocation
    block = 0
    disp_y = disp_y
    disp_x = disp_x
    gtol = 1e-5
    rtol = 1e-6
    slip_sys_file_name = input_slip_sys.txt
    C_ijkl = '279.06e3 119.6e3 119.6e3 279.06e3 119.6e3 279.06e3 79.731e3 79.731e3 79.731e3'
    nss = 48
    fill_method = symmetric9
    num_internal_var_nss = 4
    num_internal_var_scalar = 2
    intvar_read_type = slip_sys_file
    num_slip_sys_props = 4
    slip_incr_tol = 0.1
    penalty_param = 1e-5
    tan_mod_type = none
    read_prop_user_object = prop_read
    save_euler_angle = true
    burgers_length  = 0.248e-6 #mm
    disloc_line_dist = 1.5e-6  #mm
    jump_freq = 1e6
    active_enthal  = 1.2161e-18
    thermal_resist = 390.0
    temp = 323.0
    elastic_const_from_tensor = false
    maximum_substep_iteration = 10
  [../]
[]

[Postprocessors]
  [./stress_yy]
    type = ElementAverageValue
    variable = stress_yy
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./e_yy]
    type = ElementAverageValue
    variable = e_yy
    block = 'ANY_BLOCK_ID 0'
  [../]
  [./resid_y]
    type = NodalSum
    variable = resid_y
    boundary = top
  [../]
[]


[Executioner]
  type = Transient

  # Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  petsc_options_iname = -pc_hypre_type
  petsc_options_value = boomerang

  nl_rel_tol = 1e-10
  dtmax = 10.0
  end_time = 2.0
  dtmin = 1e-4
  num_steps = 10
  nl_max_its = 20

  [./TimeStepper]
   type = IterationAdaptiveDT
   reset_dt = true
   dt = 0.001
   time_t = '0.0 0.27'
   time_dt = '0.1 1e-4'
   growth_factor = 1.8
   cutback_factor = 0.5
  [../]

[]

[Outputs]
  file_base = crysp_save_euler_out
  exodus = true
  csv = true
  gnuplot = true
  [./console]
    type = Console
    perf_log = true
    output_linear = true
  [../]
[]

[Problem]
  use_legacy_uo_initialization = false
[]

[UserObjects]
  [./prop_read]
    type = ElementPropertyReadFile
    prop_file_name = 'euler_ang_file.txt'
    # Enter file data as prop#1, prop#2, .., prop#nprop
    nprop = 3
    read_type = element
  [../]
[]

