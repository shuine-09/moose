[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 10
[]

[Variables]
  [./u]
  [../]
[]

[Kernels]
  [./diff]
    type = CoefDiffusion
    variable = u
    coef = 0.1
  [../]
  [./time]
    type = TimeDerivative
    variable = u
  [../]
[]

[BCs]
  [./left]
    type = DirichletBC
    variable = u
    boundary = left
    value = 0
  [../]
  [./right]
    type = DirichletBC
    variable = u
    boundary = right
    value = 1
  [../]
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Dampers]
  [./u_damp]
    type = ConstantDamper
    variable = u
    damping = 0.9
  [../]
[]

[Outputs]
  exodus = true
[]

[Controls]
  [./damping_control]
    type = TimePeriod
    disable_objects = 'u_damp'
    start_time = 0.25
    execute_on = 'initial timestep_begin'
  [../]
[]
