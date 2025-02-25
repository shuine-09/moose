[Mesh]
  type = GeneratedMesh
  dim = 2

  xmin = 0.0
  xmax = 1.0
  nx = 10

  ymin = 0.0
  ymax = 1.0
  ny = 10
[]

[Variables]
  active = 'diffusion'

  [./diffusion]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  active = 'diff'

  [./diff]
    type = Diffusion
    variable = diffusion
  [../]
[]

[BCs]
  active = 'left right'

  [./left]
    type = DirichletBC
    variable = diffusion
    boundary = 1
    value = 3
  [../]

  [./right]
    type = DirichletBC
    variable = diffusion
    boundary = 2
    value = 1
  [../]
[]

[Dampers]
  # Use a constant damping parameter
  [./diffusion_damp]
    type = ConstantDamper
    variable = diffusion
    damping = 0.9
  [../]
[]

[Executioner]
  type = Steady

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'


[]

[Outputs]
  file_base = out
  exodus = true
[]
