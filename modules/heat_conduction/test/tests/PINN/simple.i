[Mesh]
  type = GeneratedMesh
  dim = 2
  xmax = 0.04
  ymax = 0.02
  nx = 200
  ny = 100
[]

[Variables]
  [./temp]
    initial_condition = 1000.0
  [../]
[]

[Kernels]
  [./time]
    type = HeatConductionTimeDerivative
    variable = temp
  [../]
  [./heat]
    type = HeatConduction
    variable = temp
    diffusion_coefficient = thermal_conductivity
  [../]
[]

[BCs]
  [./fix]
    type = DirichletBC
    boundary = 'left right bottom'
    variable = temp
    value = 300
  [../]
  [./flux]
    type = NeumannBC
    boundary = 'top'
    value = 0
    variable = temp
  [../]
[]

[Materials]
  [./density]
    type = GenericConstantMaterial
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
  solve_type = PJFNK
  dt = 0.01
  dtmin = 0.05
  nl_abs_tol = 1e-8
  num_steps = 100
[]


[Outputs]
  exodus = true
[]
