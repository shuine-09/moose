[Mesh]
  # type = GeneratedMesh
  # dim = 2
  # nx = 4
  # ny = 4
  # elem_type = QUAD4
  type = FileMesh
  file = crack_mesh.e
  uniform_refine = 1
  displacements = 'ux uy'
[]

[Variables]
  [./ux]
  [../]
  [./uy]
  [../]
[]

[Modules]
  [./PhaseField]
    [./Nonconserved]
      [./c]
        free_energy = E_el
        mobility = L
        kappa = kappa_op
      [../]
    [../]
  [../]
[]

[AuxVariables]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./fp_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./e_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./gss]
    order = CONSTANT
    family = MONOMIAL
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
    displacements = 'ux uy'
    use_displaced_mesh = true
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
  [../]
  [./fp_yy]
    type = RankTwoAux
    variable = fp_yy
    rank_two_tensor = fp
    index_j = 1
    index_i = 1
    execute_on = timestep_end
  [../]
  [./e_yy]
    type = RankTwoAux
    variable = e_yy
    rank_two_tensor = lage
    index_j = 1
    index_i = 1
    execute_on = timestep_end
  [../]
  [./gss]
    type = MaterialStdVectorAux
    variable = gss
    property = state_var_gss
    index = 0
    execute_on = timestep_end
  [../]
[]

[BCs]
  [./symmy]
    type = PresetBC
    variable = uy
    boundary = 1
    value = 0
  [../]
  [./symmx]
    type = PresetBC
    variable = ux
    boundary = 1
    value = 0
  [../]
  [./tdisp]
    type = FunctionPresetBC
    variable = uy
    boundary = 2
    function = tdisp
  [../]
[]

[UserObjects]
  [./slip_rate_gss]
    type = CrystalPlasticitySlipRateGSS
    variable_size = 12
    slip_sys_file_name = input_slip_sys.txt
    num_slip_sys_flowrate_props = 2
    flowprops = '1 4 0.001 0.1 5 8 0.001 0.1 9 12 0.001 0.1'
    uo_state_var_name = state_var_gss
  [../]
  [./slip_resistance_gss]
    type = CrystalPlasticitySlipResistanceGSS
    variable_size = 12
    uo_state_var_name = state_var_gss
  [../]
  [./state_var_gss]
    type = CrystalPlasticityStateVariable
    variable_size = 12
    groups = '0 4 8 12'
    group_values = '60.8 60.8 60.8'
    uo_state_var_evol_rate_comp_name = state_var_evol_rate_comp_gss
    scale_factor = 1.0
  [../]
  [./state_var_evol_rate_comp_gss]
    type = CrystalPlasticityStateVarRateComponentGSS
    variable_size = 12
    hprops = '1.0 541.5 109.8 2.5'
    uo_slip_rate_name = slip_rate_gss
    uo_state_var_name = state_var_gss
  [../]
[]

[Materials]
  [./crysp]
    type = FiniteStrainUObasedCPPF
    stol = 1e-2
    tan_mod_type = exact
    uo_slip_rates = 'slip_rate_gss'
    uo_slip_resistances = 'slip_resistance_gss'
    uo_state_vars = 'state_var_gss'
    uo_state_var_evol_rate_comps = 'state_var_evol_rate_comp_gss'
    damage_stiffness = 1e-8
    c = c
    F_name = E_el
    use_linear_fracture_energy = false
  [../]
  [./strain]
    type = ComputeFiniteStrain
    displacements = 'ux uy'
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensorCP
    C_ijkl = '1.684e5 1.214e5 1.214e5 1.684e5 1.214e5 1.684e5 0.754e5 0.754e5 0.754e5'
    fill_method = symmetric9
  [../]
  [./pfbulkmat]
    type = GenericConstantMaterial
    prop_names = 'l visco'
    prop_values = '0.02 1'
  [../]
  [./pfgc]
    type = GenericFunctionMaterial
    prop_names = 'gc_prop'
    prop_values = '1.0e-3'
  [../]
  [./define_mobility]
    type = ParsedMaterial
    material_property_names = 'gc_prop visco'
    f_name = L
    function = '1/(gc_prop * visco)'
  [../]
  [./define_kappa]
    type = ParsedMaterial
    material_property_names = 'gc_prop l'
    f_name = kappa_op
    function = 'gc_prop * l'
  [../]
[]

[Postprocessors]
  [./stress_yy]
    type = ElementAverageValue
    variable = stress_yy
  [../]
  [./fp_yy]
    type = ElementAverageValue
    variable = fp_yy
  [../]
  [./e_yy]
    type = ElementAverageValue
    variable = e_yy
  [../]
  [./gss]
    type = ElementAverageValue
    variable = gss
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
  dt = 0.01
  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
petsc_options_value = 'lu superlu_dist'

nl_rel_tol = 1e-8

l_max_its = 10
nl_max_its = 50

  dtmax = 10.0
  end_time = 1
  dtmin = 1.0e-5
  num_steps = 200
  nl_abs_step_tol = 1e-10
[]

[Outputs]
  exodus = true
  csv = true
[]
