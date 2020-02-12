[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 0.01
    ymin = 0
    ymax = 0.01
    nx = 32
    ny = 32
    elem_type = QUAD4
  []
[]

[Adaptivity]
  steps = 4
  marker = box
  max_h_level = 4
  initial_steps = 4
  stop_time = 1.0e-10
  [./Markers]
    [./box]
      bottom_left = '0.004 0.004 0'
      inside = refine
      top_right = '0.006 0.006 0'
      outside = do_nothing
      type = BoxMarker
    [../]
  [../]
[]

[ICs]
  [./ls_ic]
    type = FunctionIC
    function = ls_exact
    variable = ls
  [../]
  [velocity]
    type = VectorConstantIC
    x_value = 1e-15
    y_value = 1e-15
    variable = velocity
  []
[]

[Variables]
  [./grad_ls]
    family = LAGRANGE_VEC
  [../]
  [./curvature]
  [../]
  [./temp]
    initial_condition = 300
  [../]
[]

[AuxVariables]
  [vel_x]
  []
  [vel_y]
  []
  [density]
    family = MONOMIAL
    order = CONSTANT
  []
  [dhdT]
    family = MONOMIAL
    order = CONSTANT
  []
  [enthalpy]
    family = MONOMIAL
    order = CONSTANT
  []
  [thermal_conductivity]
    family = MONOMIAL
    order = CONSTANT
  []
  [ded_momentum_x]
    family = MONOMIAL
    order = CONSTANT
  []
  [ded_momentum_y]
    family = MONOMIAL
    order = CONSTANT
  []
  [liquid_mass_fraction]
    family = MONOMIAL
    order = CONSTANT
  []
  [solid_mass_fraction]
    family = MONOMIAL
    order = CONSTANT
  []
  [./velocity]
    family = LAGRANGE_VEC
    order = FIRST
  [../]
  [./ls]
    order = FIRST
  [../]
[]

[Functions]
  [./ls_exact]
    type = LevelSetOlssonPlane
    epsilon = 0.0002
  [../]
[]

[AuxKernels]
  [vel_x]
    type = VectorVariableComponentAux
    variable = vel_x
    vector_variable = velocity
    component = 'x'
  []
  [vel_y]
    type = VectorVariableComponentAux
    variable = vel_y
    vector_variable = velocity
    component = 'y'
  []
  [density]
    type = MaterialRealAux
    variable = density
    property = rho
    execute_on = timestep_end
  []
  [dhdT]
    type = MaterialRealAux
    variable = dhdT
    property = dhdT
    execute_on = timestep_end
  []
  [enthalpy]
    type = MaterialRealAux
    variable = enthalpy
    property = enthalpy
    execute_on = timestep_end
  []
  [thermal_conductivity]
    type = MaterialRealAux
    variable = thermal_conductivity
    property = thermal_conductivity
    execute_on = timestep_end
  []
  [ded_momentum_x]
    type = MaterialRealVectorValueAux
    variable = ded_momentum_x
    property = ded_momentum
    component = 0
    execute_on = timestep_end
  []
  [ded_momentum_y]
    type = MaterialRealVectorValueAux
    variable = ded_momentum_y
    property = ded_momentum
    component = 1
    execute_on = timestep_end
  []
  [liquid_mass_fraction]
    type = MaterialRealAux
    variable = liquid_mass_fraction
    property = liquid_mass_fraction
    execute_on = timestep_end
  []
  [solid_mass_fraction]
    type = MaterialRealAux
    variable = solid_mass_fraction
    property = solid_mass_fraction
    execute_on = timestep_end
  []
[]

[Kernels]
  [./grad_ls]
    type = LevelSetGradientRegulization
    level_set = ls
    variable = grad_ls
  [../]
  [./curvature]
    type = LevelSetCurvatureRegulization
    grad_level_set = grad_ls
    variable = curvature
    epsilon = 2e-7
  [../]
  [./heat_time]
    type = ADHeatConductionTimeDerivative
    specific_heat = dhdT
    density_name = rho
    variable = temp
  [../]
  [./heat_conv]
    type = HeatConvection
    variable = temp
    velocity = velocity
    density = 'rho'
    enthalpy = 'enthalpy'
  [../]
  [./heat_cond]
    type = ADHeatConduction
    thermal_conductivity = thermal_conductivity
    variable = temp
  [../]
  [./heat_source]
    type = DEDHeatSource
    level_set = ls
    grad_level_set = grad_ls
    variable = temp
    rho_m = rho_ixture
    rho = rho
    h_m = enthalpy_mixture
    h = enthalpy
    rho_g = 1.184
    c_g = 1000
    laser_power = 300
    effective_beam_radius = 5e-5
    absorption_coefficient = 0.27
    heat_transfer_coefficient = 100
    StefanBoltzmann_constant = 5.67e-8
    material_emissivity = 0.59
    ambient_temperature = 300
    laser_center_x = '0.005+12.7*1e-3*t'
    laser_center_y = '0.005'
  [../]
[]

[Materials]
  [ded_mat]
    type = DEDMaterial
    level_set = ls
    grad_level_set = grad_ls
    temperature = temp
    curvature = curvature
    velocity = velocity
    rho_g = 1.184
    rho_s = 8000
    rho_l = 8000
    mu_g = 1.81e-5
    mu_l = 0.1
    mu_s = 10
    c_g = 1000
    c_s = 500
    c_l = 500
    k_g = 26.24e-3
    k_s = 31.8724
    k_l = 209.3
    solidus_temperature = 1648
    liquidus_temperature = 1673
    latent_heat = 2.5e5
    capillary_coefficient = 1.169
    thermalcapillary_coefficient = -3e-4#-4.3e-4
    laser_center_x = '0.005+12.7*1e-3*t'
    laser_center_y = '0.005'
    permeability_constant = 1e-6
    gravity = '0 -9.8 0'
    thermal_expansion = 1.45e-4
  []
[]

[BCs]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
    solve_type = 'NEWTON'
  [../]
[]

[Executioner]
  #automatic_scaling = true
  type = Transient
  solve_type = NEWTON
  start_time = 0
  dt = 1e-4
  end_time = 0.025
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu superlu_dist'
  #line_search = 'cp'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-6
  nl_max_its = 10
  l_tol = 1e-6
  l_max_its = 20
  #dtmin = 1e-3
  dtmax = 1
  line_search = 'none'
[]

[Outputs]
  exodus = true
[]
