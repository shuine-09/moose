[Mesh]
  type = GeneratedMesh
  dim = 2
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1
  nx = 20
  ny = 20
  uniform_refine = 3 #1/64
[]

[Problem]
  coord_type = RZ
[]

[AuxVariables]
  [./velocity]
    family = LAGRANGE_VEC
  [../]
[]

[Variables]
  [./phi]
  [../]
[]

[Functions]
  [./phi_exact]
    type = LevelSetOlssonPlane
    epsilon = 0.01
    point = '0.5 0.5 0'
    normal = '0 1 0'
  [../]
  [./velocity_func]
    type = ParsedVectorFunction
    value_x = 'y'
    value_y = '10*(1-x)^2'
  [../]
[]

[ICs]
  [./phi_ic]
    type = FunctionIC
    function = phi_exact
    variable = phi
  [../]
  [./vel_ic]
    type = VectorFunctionIC
    variable = velocity
    function = velocity_func
  []
[]

[Kernels]
  [./time]
    type = TimeDerivative
    variable = phi
  [../]

  [./advection]
    type = LevelSetAdvection
    velocity = velocity
    variable = phi
  [../]
[]

[Postprocessors]

  [./cfl]
    type = LevelSetCFLCondition
    velocity = velocity
    execute_on = 'initial'
  [../]

[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  start_time = 0
  end_time = 1
  nl_rel_tol = 1e-12
  scheme = crank-nicolson
  petsc_options_iname = '-pc_type -pc_sub_type'
  petsc_options_value = 'asm      ilu'
  [./TimeStepper]
    type = PostprocessorDT
    postprocessor = cfl
    scale = 1
  [../]

[]

[MultiApps]
  [./reinit]
    type = LevelSetReinitializationMultiApp
    input_files = 're.i'
    execute_on = 'timestep_end'
  [../]
[]

[Transfers]
  [./to_sub]
    type = MultiAppCopyTransfer
    variable = phi
    source_variable = phi
    direction = to_multiapp
    multi_app = reinit
    execute_on = 'timestep_end'
  [../]

  [./to_sub_init]
    type = MultiAppCopyTransfer
    variable = phi_0
    source_variable = phi
    direction = to_multiapp
    multi_app = reinit
    execute_on = 'timestep_end'
  [../]

  [./from_sub]
    type = MultiAppCopyTransfer
    variable = phi
    source_variable = phi
    direction = from_multiapp
    multi_app = reinit
    execute_on = timestep_end
  [../]
[]

[Outputs]
  exodus = true
  csv = true
[]
