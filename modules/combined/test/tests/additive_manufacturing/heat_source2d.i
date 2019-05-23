[Mesh]
  type = MeshGeneratorMesh
  uniform_refine = 0
[]

[MeshGenerators]
  [./mesh]
    type = GeneratedMeshGenerator
    nx = 100
    ny = 20
    xmax = 100
    ymax = 20
    dim = 2
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

[BCs]
  [./convective_substrate]
    type = ConvectiveFluxFunction # Convective flux, e.g. q'' = h*(Tw - Tf)
    variable = temp
    boundary = 'left right top bottom'
    coefficient = 1
    T_infinity = 300
  [../]
[]

[Materials]
  [./heat]
    type = ADHeatConductionMaterial
    specific_heat = 603
    thermal_conductivity = 10e-2
  [../]
  [./volumetric_heat]
    type = DoubleEllipsoidHeatSource
    a = 16
    b = 8
    c = 1.5
    power = 800
    efficienty = 0.45
    factor = 2
    velocity = 10
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

  dt = 0.01
  end_time = 10
[]

[Outputs]
  exodus = true
[]
