#This input does not add time derivative kernel for phase field equation
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 800
  ny = 799
  xmin = 0
  xmax = 4
  ymin = 0
  ymax = 4
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
  out_of_plane_strain = strain_zz
[]

[XFEM]
  qrule = volfrac
  output_cut_plane = true
[]

[Variables]
  [./strain_zz]
  [../]
[]

[Modules]
  [./TensorMechanics]
    [./Master]
      [./mech]
        add_variables = true
        strain = SMALL
        additional_generate_output = 'stress_yy stress_zz stress_xx stress_xy strain_xx strain_yy strain_zz strain_xy'
        planar_formulation = WEAK_PLANE_STRESS
      [../]
    [../]
  [../]
[]

[UserObjects]
  [./line_seg_cut_uo]
    type = LineSegmentCutUserObject
    cut_data = '1.8 2 2.2 2'
  [../]
[]


[BCs]
  [./yfix]
    type = PresetBC
    variable = disp_y
    boundary = 'left right top bottom'
    value = 0
  [../]
  [./xfix]
    type = PresetBC
    variable = disp_x
    boundary = 'left right top bottom'
    value = 0
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    C_ijkl = '1 0.2'
    fill_method = symmetric_isotropic_E_nu
  [../]
  [./elastic]
    type = ComputeLinearElasticStress
  [../]
[]

[DiracKernels]
  [./pressure_x]
    type = XFEMPressure
    variable = disp_x
    component = 0
    function = 1.0e-3
  [../]

  [./pressure_y]
    type = XFEMPressure
    variable = disp_y
    component = 1
    function = 1.0e-3
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

  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu       superlu_dist'

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-6
  l_max_its = 10
  nl_max_its = 50

  dt = 1e-4
  num_steps = 3
[]

[Outputs]
  exodus = true
  csv = true
[]
