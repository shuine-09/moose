[Problem]
  XFEM_cuts = '0.55   0.399   0.55   0.601   0.0000e+00   1.0000e+00'
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 5
  ny = 5
  xmin = 0.0;
  xmax = 1.0;
  ymin = 0.0;
  ymax = 1.0;
  nz = 0
  zmax = 0
  elem_type = QUAD4
[]

[AuxVariables]
  [./xfem_volfrac]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[CutPlanes]
  order = CONSTANT
  family = MONOMIAL
[]

[AuxKernels]
  [./xfem_volfrac]
    type = XFEMVolFracAux
    variable = xfem_volfrac
    execute_on = timestep_begin
  [../]
[]

[Variables]
  [./u]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Functions]
  [./force]
    type = ParsedFunction
    value = t
  [../]
[]

[Kernels]
  [./diff]
    type = Diffusion
    variable = u
  [../]
  [./force]
    type = UserForcingFunction
    variable = u
    function = force
  [../]
[]

[BCs]
  [./right]
    type = DirichletBC
    variable = u
    boundary = right
    value = 1
  [../]
  [./left]
    type = DirichletBC
    variable = u
    boundary = left
    value = 0
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 4
  dt = 1

  # Preconditioned JFNK (default)
  solve_type = 'PJFNK'

[]

[Adaptivity]
  steps = 1
  marker = box
  max_h_level = 1
  initial_steps = 1
  [./Markers]
    [./box]
      bottom_left = '0.2 0.2 0'
      inside = refine
      top_right = '0.8 0.8 0'
      outside = do_nothing
      type = BoxMarker
    [../]
  [../]
[]

[Outputs]
  exodus = true
  output_initial = true
  print_perf_log = true
[]
