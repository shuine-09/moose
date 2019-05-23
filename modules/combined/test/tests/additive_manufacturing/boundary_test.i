[Mesh]
  type = FileMesh
  file = boundary_test.e
  displacements = 'disp_x disp_y'
[]

[GlobalParams]
  block = '1 2'
[]

[Variables]
  [./temp]
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
  [./null]
    type = NullKernel
    variable = temp
    block = 2
  [../]
[]

[BCs]
  [./bottom_temp]
    type = DirichletBC
    variable = temp
    boundary = '101'
    value = 10.0
  [../]
  [./interior]
    type = PenaltyDirichletBC
    variable = temp
    boundary = 1
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
  type = Transient
  solve_type = 'NEWTON'

  nl_rel_tol = 1e-14
  l_tol = 1e-3

  l_max_its = 100

  dt = 1.0
  end_time = 1.0
[]

[Outputs]
  exodus = true
[]
