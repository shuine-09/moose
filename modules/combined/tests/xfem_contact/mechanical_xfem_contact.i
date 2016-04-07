# This is a mechanical constraint (contact formulation) problem that pushes a square against a cracked element
# Switch penalty from 0 to 1e6 to compare contact vs non-contact

[Mesh]
  file = elements.msh
  displacements = 'disp_x disp_y'
[]


[XFEM]
  cut_data = '2.1 0.5 1.1 0.5 0 0'
  qrule = volfrac
  output_cut_plane = true
[]



[Variables]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]

  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]

  [] # Variables

[AuxVariables]

  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]

[] # AuxVariables

[SolidMechanics]
  [./solid]
    disp_x = disp_x
    disp_y = disp_y
  [../]
[]

[Contact]
  [./dummy_name]
    master = 3
    slave = 4
    disp_x = disp_x
    disp_y = disp_y
    penalty = 1e6 # contact enforced
    #penalty = 0e6 # no contact
    formulation = kinematic
    system = constraint
  [../]
[]


[AuxKernels]
  [./stress_xx]
    type = MaterialTensorAux
    tensor = stress
    variable = stress_xx
    index = 0
  [../]
[] # AuxKernels

[BCs]
  [./left_x]
    type = DirichletBC
    variable = disp_x
    boundary = 1
    #value = 0.05
    value = 0.2
  [../]

  [./left_y]
    type = DirichletBC
    variable = disp_y
    boundary = 1
    value = 0.0
  [../]



  [./right_x]
    type = DirichletBC
    variable = disp_x
    boundary = 2
    value = 0.0
  [../]

  [./right_y]
    type = DirichletBC
    variable = disp_y
    boundary = 2
    value = 0.0
  [../]



[] # BCs

[Materials]

  [./stiffStuff1]
    type = Elastic
    block = 1

    disp_x = disp_x
    disp_y = disp_y


    youngs_modulus = 1e6
    poissons_ratio = 0.3
  [../]
  [./stiffStuff2]
    type = Elastic
    block = 2

    disp_x = disp_x
    disp_y = disp_y


    youngs_modulus = 1e6
    poissons_ratio = 0.3
  [../]
[] # Materials

[Executioner]
  type = Transient

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'



  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre    boomeramg      101'


  line_search = 'none'


  nl_abs_tol = 1e-8

  l_max_its = 100
  nl_max_its = 10
  dt = 1.0
  num_steps = 1
[] # Executioner

[Outputs]
  file_base = mechanical_constraint_out_comparison1
  [./exodus]
    type = Exodus
     #exodus = true
    elemental_as_nodal = true
  [../]
[] # Outputs
