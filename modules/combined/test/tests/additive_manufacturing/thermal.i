[Mesh]
  type = GeneratedMesh
  nx = 10
  ny = 10
  dim = 2
  elem_type = QUAD4
  displacements = 'disp_x disp_y'
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
  [../]
  [./heat_source]
    type = ADMatHeatSource
    material_property = volumetric_heat
    variable = temp
    scalar = 10
    use_displaced_mesh = true
  [../]
[]

[BCs]
  [./bottom_temp]
    type = DirichletBC
    variable = temp
    boundary = 1
    value = 10.0
  [../]
[]

[Materials]
  [./heat]
    type = HeatConductionMaterial
    specific_heat = 1.0
    thermal_conductivity = 1.0
  [../]
  [./volumetric_heat]
    type = GenericConstantMaterial
    prop_names = 'volumetric_heat'
    prop_values = '3.8e1'
  [../]
[]

[MultiApps]
  [disp]
    type = TransientMultiApp
    positions = '0.0 0.0 0.0'
    input_files = disp.i
    # sub_cycling = true
    # steady_state_tol = 1e-6
    # detect_steady_state = true
    execute_on = 'timestep_end'
  []
[]

[Transfers]
  [to_temp]
    type = MultiAppMeshFunctionTransfer
    direction = to_multiapp
    execute_on = 'timestep_end'
    multi_app = disp
    source_variable = temp
    variable = temp_aux
  []
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
