[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Mesh]
  type = FileMesh
  file = crack_mesh.e
  uniform_refine = 2
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./c]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxVariables]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./dcx]
    order = FIRST
    family = MONOMIAL
  [../]
  [./dcy]
    order = FIRST
    family = MONOMIAL
  [../]
[]

[Kernels]
  [./pfbulk]
    type = AllenCahn
    variable = c
    mob_name = M
    f_name = E_el
  [../]
  [./TensorMechanics]
    displacements = 'disp_x disp_y'
  [../]
  # [./solid_x]
  #   type = PhaseFieldFractureMechanicsOffDiag
  #   variable = disp_x
  #   component = 0
  #   c = c
  # [../]
  # [./solid_y]
  #   type = PhaseFieldFractureMechanicsOffDiag
  #   variable = disp_y
  #   component = 1
  #   c = c
  # [../]
  # [./dcdt]
  #   type = CoefTimeDerivative
  #   variable = c
  #   Coefficient = 1e-5
  # [../]
  [./dcdt]
    type = TimeDerivative
    variable = c
    # Coefficient = 1e-6
  [../]
  [./InterfacialE]
    type = AnisotropicGradEnergy
    variable = c
    mob_name = M
    kappa_name = Wsq_aniso
    gradient_component_names = 'dcx dcy'
    # gradmag_threshold = 1e-10  #threshold
  [../]
  [./Rc_aniso_doublewell]
    type = AnisotropicDoubleWellEnergy
    variable = c
    mob_name = M
    fbulk_name = f_aniso_m4
    gradient_component_names = 'dcx dcy'
    # gradmag_threshold = 1e-10
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
  [./get_dcx]
    type = VariableGradientComponent
    variable = dcx
    gradient_variable = c
    component = x
  [../]
  [./get_dcy]
    type = VariableGradientComponent
    variable = dcy
    gradient_variable = c
    component = y
  [../]
[]

[BCs]
  [./ydisp]
    type = FunctionPresetBC
    variable = disp_y
    boundary = 2
    function = 't'
  [../]
  [./xfixleft]
    type = PresetBC
    variable = disp_x
    boundary = 3
    value = 0
  [../]
  [./yfix]
    type = PresetBC
    variable = disp_y
    boundary = '1'
    value = 0
  [../]
  # [./xfix]
  #   type = PresetBC
  #   variable = disp_x
  #   boundary = 1
  #   value = 0
  # [../]
[]

[Materials]
  [./pfbulkmat]
    type = GenericConstantMaterial
    prop_names = 'gc_prop gc l visco eps_bar delta M'
    prop_values = '0    1e-3 0.02 1.0 4e-2 0.4 1e6'
  [../]
  [./mat]
    type = DerivativeParsedMaterial
    f_name = Wsq_aniso
    material_property_names = 'gc delta l cos2theta(dcx,dcy) sin2theta(dcx,dcy) cos2theta1(dcx,dcy) sin2theta1(dcx,dcy)'
    args = 'dcx dcy'
    #function = 'if(dcx^2 + dcy^2 > 1e-4, gc*l*(1 + delta * (cos2theta1*cos2theta0 + sin2theta1*sin2theta0))^2, gc*l)'  #use threshold
    function = 'gc*l*(1 + delta * (cos2theta*cos2theta0 + sin2theta*sin2theta0))^2'
    constant_names       = 'sin2theta0  cos2theta0'
    constant_expressions = '0.8660254 0.5'
    derivative_order = 2
  [../]
  [./f_aniso_m4]
    type = DerivativeParsedMaterial
    f_name = f_aniso_m4
    # eps4 is the anisotropy strength
    material_property_names = 'gc delta l cos2theta(dcx,dcy) sin2theta(dcx,dcy) cos2theta1(dcx,dcy) sin2theta1(dcx,dcy)'
    # constant_names = 'e'
    # constant_expressions = '0.01'
    args = 'c dcx dcy'
    #function = 'if(dcx^2 + dcy^2 > 1e-4, gc/l*c*c*0.5*(1 + delta * (cos2theta1*cos2theta0 + sin2theta1*sin2theta0))^2, gc/l*c*c*0.5)' ##use threshold
    function = 'gc/l*c*c*0.5*(1 + delta * (cos2theta*cos2theta0 + sin2theta*sin2theta0))^2'
    constant_names       = 'sin2theta0  cos2theta0'
    constant_expressions = '0.8660254 0.5'
    derivative_order = 2
    # outputs = exodus
  [../]
  [./sin2theta1]
    type = DerivativeParsedMaterial
    f_name = sin2theta1
    args = 'dcx dcy'
    constant_names = 'e'
    constant_expressions = '1e-2'
    function = '2.0*dcx*dcy/(dcx^2 + dcy^2)'
    derivative_order = 2
  [../]
  [./cos2theta1]
    type = DerivativeParsedMaterial
    f_name = cos2theta1
    args = 'dcx dcy'
    constant_names = 'e'
    constant_expressions = '1e-2'
    function = '(dcx^2-dcy^2)/(dcx^2 + dcy^2)'
  [../]
  [./sin2theta]
    type = DerivativeParsedMaterial
    f_name = sin2theta
    args = 'dcx dcy'
    constant_names = 'e'
    constant_expressions = '1e-2'
    function = '2.0*dcx*dcy/(dcx^2 + dcy^2+e^2)'
    derivative_order = 2
  [../]
  [./cos2theta]
    type = DerivativeParsedMaterial
    f_name = cos2theta
    args = 'dcx dcy'
    constant_names = 'e'
    constant_expressions = '1e-2'
    function = '(dcx^2-dcy^2)/(dcx^2 + dcy^2+e^2)'
  [../]
  [./elastic]
    type = ComputeLinearElasticPFFractureStress
    #type = ComputeIsotropicLinearElasticPFFractureStress
    kdamage = 1e-6
    c = c
    F_name = E_el
    use_current_history_variable = true
    use_current_damage = false
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    C_ijkl = '120.0 80.0'
    fill_method = symmetric_isotropic
  [../]
  [./strain]
    type = ComputeSmallStrain
  [../]
  [./Frac_density]
    type = ParsedMaterial
    f_name = fFrac
    material_property_names = 'gc_prop l'
    args = 'c dcx dcy'
    function = '0.5*gc_prop*l*(dcx^2+dcy^2) + 0.5*gc_prop*c*c/l'
  [../]
  [./interfacial_density]
    type = ParsedMaterial
    f_name = finterf
    material_property_names = 'gc_prop l'
    args = 'dcx dcy'
    function = '0.5*gc_prop*l*(dcx^2+dcy^2)'
  [../]
  [./tot_free_density]
    type = ParsedMaterial
    f_name = ftot
    material_property_names = 'gc_prop l E_el'
    args = 'c dcx dcy'
    function = '0.5*gc_prop*l*(dcx^2+dcy^2) + E_el'
  [../]
[]

[Preconditioning]
  active = 'SMP'
  [./SMP]
    type = SMP
    full = true
  [../]
  [./FDP]
    type = FDP
    full = true
  [../]
[]
[Postprocessors]
  [./Total_frac_energy]
    type = ElementIntegralMaterialProperty
    mat_prop = fFrac
    outputs = csv
  [../]
  [./Interfacial_energy]
    type = ElementIntegralMaterialProperty
    mat_prop = finterf
    outputs = csv
  [../]
  [./total_free_energy]
    type = ElementIntegralMaterialProperty
    mat_prop = ftot
    outputs = csv
  [../]
  [./disp_y_top]
    type = SideAverageValue
    variable = disp_y
    boundary = 2
  [../]
  [./stressyy]
    type = ElementAverageValue
    variable = stress_yy
  [../]
[]

[Executioner]
  type = Transient

  solve_type = NEWTON

  scheme = bdf2

  petsc_options_iname = '=ksp_type -pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'preonly lu superlu_dist'

  line_search = 'none'

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-8
  l_tol = 1e-4
  l_max_its = 30
  nl_max_its = 15
  # [./Adaptivity]
  #   initial_adaptivity = 0
  #   refine_fraction = 0.95
  #   coarsen_fraction = 0.01
  #   max_h_level = 0
  #   weight_names = c
  #   weight_values = 1
  # [../]
  dt = 1e-5
  dtmin = 1e-10
  start_time = 0.0
  end_time = 0.008
[]

[Outputs]
  file_base = R30
  exodus = true
  csv = true
  gnuplot = true
  interval = 1
[]
