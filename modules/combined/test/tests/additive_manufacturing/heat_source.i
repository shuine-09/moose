[Mesh]
  type = MeshGeneratorMesh
  uniform_refine = 1
[]

[MeshGenerators]
  [./mesh]
    type = GeneratedMeshGenerator
    nx = 3
    ny = 3
    nz = 38
    xmin = -1.5
    xmax = 1.5
    ymin = -1.5
    ymax = 1.5
    zmin = 0
    zmax = 38
    dim = 3
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
[]

[Kernels]
  [./time]
    type = ADHeatConductionTimeDerivative
    variable = temp
  [../]
  [./heat_conduct]
    type = ADHeatConduction
    variable = temp
    use_displaced_mesh = true
  [../]
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

[Materials]
  [./heat]
    type = ADHeatConductionMaterial
    specific_heat = 603
    thermal_conductivity = 10e-3
  [../]
  [./volumetric_heat]
    type = DoubleEllipsoidHeatSource
    a = 1.5
    b = 1.5
    c = 1.5
    power = 425
    efficienty = 0.45
    factor = 1
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

  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-12
  l_tol = 1e-3

  l_max_its = 100

  dt = 0.0885827
  end_time = 4.0
[]

[Outputs]
  exodus = true
[]
