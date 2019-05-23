[Mesh]
  type = FileMesh
  file = boundary_move.e
  displacements = 'disp_x disp_y'
[]

[MeshModifiers]
  [./translate]
    type = Transform
    transform = TRANSLATE
    vector_value = '5 0 0'
  [../]
[]

[Variables]
  [./temp]
    initial_condition = 300
  [../]
[]

[AuxVariables]
  [./disp_x]
    initial_condition = 0
  [../]
  [./disp_y]
    initial_condition = 0
  [../]
[]

[Kernels]
  [./time]
    type = TimeDerivative
    variable = temp
  [../]
  [./heat_conduct]
    type = ADHeatConduction
    variable = temp
    use_displaced_mesh = true
    block = 1
  [../]
  [./heat_source]
    type = ADMatHeatSource
    material_property = volumetric_heat
    variable = temp
    scalar = 10
    use_displaced_mesh = true
    block = 1
  [../]
[]

[BCs]
  [./bottom_temp]
    type = DirichletBC
    variable = temp
    boundary = '101'
    value = 300.0
  [../]
  [./interior]
    type = PenaltyDirichletBC
    variable = temp
    boundary = 102
    value = 100
    penalty = 1e10
  [../]
[]

[Materials]
  [./heat]
    type = ADHeatConductionMaterial
    specific_heat = 1.0
    thermal_conductivity = 1.0
  [../]
  [./volumetric_heat]
    type = ADGenericConstantMaterial
    prop_names = 'volumetric_heat'
    prop_values = '3.8e1'
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  automatic_scaling = true
  type = Transient
  solve_type = 'NEWTON'

  nl_rel_tol = 1e-16
  l_tol = 1e-3

  l_max_its = 100

  dt = 1.0
  end_time = 10.0
[]

[Outputs]
  exodus = true
[]
