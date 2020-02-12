[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 2
    ymin = 0
    ymax = 2
    nx = 32
    ny = 32
    elem_type = QUAD9
  []
[]

[AuxVariables]
  [vel_x]
  []
  [vel_y]
  []
[]

[ICs]
  [./phi_ic]
    type = FunctionIC
    function = phi_exact
    variable = phi
  [../]
[]

[Variables]
  [./phi]
    order = SECOND
  [../]
  [temp]
    order = SECOND
  []
[]

[Problem]
  type = LevelSetProblem
[]

[Functions]
  [./phi_exact]
    type = LevelSetOlssonPlane
    epsilon = 0.04
  [../]
[]

[Kernels]
  [./time]
    type = TimeDerivative
    variable = phi
  [../]
  [./advection]
    type = LevelSetAdvection
    velocity_x = vel_x
    velocity_y = vel_y
    variable = phi
  [../]
  [./mass_addition]
    type = LevelSetDEDMass
    variable = phi
  [../]
  [./heat_time]
    type = ADHeatConductionTimeDerivative
    specific_heat = dhdT
    density_name = rho
    variable = temp
  [../]
  # [./heat_conv]
  #   type = HeatConvection
  #   variable = temp
  #   velocity =
  # [../]
  [./heat_cond]
    type = ADHeatConduction
    thermal_conductivity = thermal_conductivity
    variable = temp
  [../]
  [./heat_source]
    type = DEDHeatSource
    variable = temp
    level_set = phi
    rho_m = rho_ixture
    rho = rho
    h_m = enthalpy_mixture
    h = enthalpy
    rho_g = 1
    c_g = 1
    laser_power = 1
    effective_beam_radius = 0.1
    absorption_coefficient = 1
    heat_transfer_coefficient = 1
    StefanBoltzmann_constant = 1
    material_emissivity = 1
    ambient_temperature = 0
    laser_center_x = '10*t'
    laser_center_y = '1'
  [../]
[]

[Materials]
  [ded_mat]
    type = DEDMaterial
    level_set = phi
    temperature = temp
    rho_g = 1
    rho_s = 1
    rho_l = 1
    mu_g = 1
    mu_s = 1
    mu_l = 1
    c_g = 1
    c_s = 1
    c_l = 1
    k_g = 1
    k_s = 1
    k_l = 1
    solidus_temperature = 0
    liquidus_temperature = 1
    latent_heat = 0
    laser_center_x = '10*t'
    laser_center_y = '1'
  []
[]

[BCs]
  [./temp]
    type = ADDirichletBC
    variable = temp
    value = 0
    boundary = 'bottom right left top'
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
    solve_type = 'NEWTON'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  start_time = 0
  dt = 0.001
  end_time = 1000
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu mumps'
  line_search = 'none'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-10
  nl_max_its = 15
  l_tol = 1e-6
  l_max_its = 20
  dtmax = 1
[]

[Outputs]
  exodus = true
[]
