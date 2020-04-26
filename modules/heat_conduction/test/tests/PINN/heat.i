[Mesh]
  type = GeneratedMesh
  dim = 2
  xmax = 0.04
  ymax = 0.02
  nx = 200
  ny = 100
[]

[GlobalParams]
  heat_convection_coef = 100
  ambient_temperature = 300
  emissivity = 0.4
  stefan_boltzmann = 5.67e-8
  substrate_thickness = 0.002
  absorptivity = 0.4
  laser_power = 500
  laser_radius = 0.00065
  laser_velocity = 0.0066
[]

[Variables]
  [./temp]
    initial_condition = 300.0
  [../]
[]

[Kernels]
  [./time]
    type = ADHeatConductionTimeDerivative
    variable = temp
  [../]
  [./heat]
    type = ADHeatConduction
    variable = temp
    thermal_conductivity = thermal_conductivity
  [../]
  [./heatsource]
    type = ADMatHeatSource
    variable = temp
  [../]
[]

[BCs]
  [./flux0]
    type = ADConvectiveHeatFluxBC
    boundary = 'left right bottom'
    variable = temp
    absorptivity = 0
  [../]
  [./flux]
    type = ADConvectiveHeatFluxBC
    boundary = 'top'
    variable = temp
  [../]
[]

[Materials]
  [./density]
    type = ADGenericConstantMaterial
    prop_names = 'density  thermal_conductivity specific_heat'
    prop_values = '8000 200 500'
  [../]
[]

[Preconditioning]
  [./full]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  dt = 0.01
  dtmin = 0.05
  nl_abs_tol = 1e-10
  num_steps = 100
[]


[Outputs]
  exodus = true
[]
