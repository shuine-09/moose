[Mesh]
  type = MeshGeneratorMesh
[]

[MeshGenerators]
  [./mesh]
    type = GeneratedMeshGenerator
    nx = 2
    ny = 2
    nz = 2
    xmax = 10
    ymax = 4
    zmax = 4
    dim = 3
  [../]
  [./add_bnd]
    type = SideSetsFromAllElementFaces
    input = mesh
    block = '0'
    new_boundary = 'moving_boundary'
  [../]
[]

[Variables]
  [./temp]
    initial_condition = 300
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
    block = 0
  [../]
  [./heat_source]
    type = ADMatHeatSource
    material_property = volumetric_heat
    variable = temp
    scalar = 10
    use_displaced_mesh = true
    block = 0
  [../]
[]

[BCs]
  [./left_temp]
    type = DirichletBC
    variable = temp
    boundary = 'left'
    value = 300.0
  [../]
  [./interior]
    type = PenaltyDirichletBC
    variable = temp
    boundary = 'moving_boundary'
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
    prop_values = '0'
  [../]
  [./density]
    type = ADDensity
    density = 1
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

  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-12
  l_tol = 1e-3

  l_max_its = 100

  dt = 1.0
  end_time = 5
[]

[Outputs]
  exodus = true
[]
