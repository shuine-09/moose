[Mesh]
  type = MeshGeneratorMesh
  uniform_refine = 0
[]

[MeshGenerators]
  [./mesh]
    type = FileMeshGenerator
    file = hollow.e
  [../]
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
  [./disp_z]
    initial_condition = 0
  [../]
[]

[Kernels]
  [./time]
    type = ADHeatConductionTimeDerivative
    variable = temp
  [../]
  # [./heat_conduct]
  #   type = ADHeatConduction
  #   variable = temp
  #   use_displaced_mesh = true
  # [../]
  [./heat_source]
    type = ADMatHeatSource
    material_property = volumetric_heat
    variable = temp
    scalar = 1
    use_displaced_mesh = true
  [../]
[]

# [BCs]
#   [./bottom_temp]
#     type = DirichletBC
#     variable = temp
#     boundary = 1
#     value = 10.0
#   [../]
# []

[MultiApps]
  [thermo_mech]
    type = TransientMultiApp
    positions = '0.0 0.0 0.0'
    input_files = sub.i
    # sub_cycling = true
    # steady_state_tol = 1e-6
    # detect_steady_state = true
    execute_on = 'timestep_end'
  []
[]

[Transfers]
  [to_mech]
    type = MultiAppCopyTransfer
    direction = to_multiapp
    execute_on = 'timestep_end'
    multi_app = thermo_mech
    source_variable = temp
    variable = temp_aux
  []
  [to_disp_x]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    execute_on = 'timestep_end'
    multi_app = thermo_mech
    source_variable = disp_x
    variable = disp_x
  []
  [to_disp_y]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    execute_on = 'timestep_end'
    multi_app = thermo_mech
    source_variable = disp_y
    variable = disp_y
  []
  [to_disp_z]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    execute_on = 'timestep_end'
    multi_app = thermo_mech
    source_variable = disp_z
    variable = disp_z
  []
 [to_temp]
   type = MultiAppCopyTransfer
   direction = from_multiapp
   execute_on = 'timestep_end'
   multi_app = thermo_mech
   source_variable = temp
   variable = temp
 []
[]

[Materials]
  [./heat]
    type = ADHeatConductionMaterial
    specific_heat = 603
    thermal_conductivity = 10e-2
  [../]
  [./volumetric_heat]
    type = DoubleEllipsoidHeatSource
    a = 5
    b = 5
    c = 3
    power = 1000
    efficienty = 0.5
    factor = 2
    velocity = -8.47
  [../]
  [./density]
    type = ADDensity
    density = 4.43e-6
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

# [Adaptivity]
#   [./Indicators]
#     [./jump]
#       type = ValueJumpIndicator
#       variable = temp
#       outputs = exodus
#     [../]
#   [../]
#   [./Markers]
#     [./error]
#       type = ErrorFractionMarker
#       indicator = jump
#       coarsen = 0.1
#       refine = 0.95
#     [../]
#   [../]
#   marker = error
#   max_h_level = 2
#   cycles_per_step = 2
# []

[Executioner]
  type = Transient
  solve_type = 'NEWTON'

  petsc_options_iname = '-ksp_type -pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'preonly lu       superlu_dist'

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6
  l_tol = 1e-3

  l_max_its = 100

  dt = 0.01
  end_time = 10
[]

[Outputs]
  exodus = true
[]
