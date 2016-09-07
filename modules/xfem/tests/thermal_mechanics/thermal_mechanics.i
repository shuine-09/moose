[GlobalParams]
  order = FIRST
  family = LAGRANGE
  disp_x = disp_x
  disp_y = disp_y
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 11
  xmin = 0.0
  xmax = 1.0
  ymin = 0.0
  ymax = 1.0
  elem_type = QUAD4
  displacements = 'disp_x disp_y'
[]

[XFEM]
  cut_data = '0.0 0.5 0.5 0.5 0 0'
  qrule = volfrac
  output_cut_plane = true
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./temp]
    initial_condition = 300.0
  [../]
[]

[AuxVariables]
  [./stress_xx]      # stress aux variables are defined for output; this is a way to get integration point variables to the output file
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[SolidMechanics]
  [./solid]
    temp = temp
  [../]
[]

[Kernels]
  [./heat]         # gradient term in heat conduction equation
    type = HeatConduction
    variable = temp
#    use_displaced_mesh = true
  [../]
[]

[AuxKernels]
  [./stress_xx]               # computes stress components for output
    type = MaterialTensorAux
    tensor = stress
    variable = stress_xx
    index = 0
    execute_on = timestep_end     # for efficiency, only compute at the end of a timestep
  [../]
  [./stress_yy]
    type = MaterialTensorAux
    tensor = stress
    variable = stress_yy
    index = 1
    execute_on = timestep_end
  [../]
[]

[Functions]
  [./pull]
    type = PiecewiseLinear
    x='0  2 4'
    y='0  0.02 0.04'
  [../]
[]

[BCs]
  [./bottomx]
    type = PresetBC
    boundary = bottom
    variable = disp_x
    value = 0.0
  [../]
  [./bottomy]
    type = PresetBC
    boundary = bottom
    variable = disp_y
    value = 0.0
  [../]
  [./topx]
    type = PresetBC
    boundary = top
    variable = disp_y
    value = 0.0
  [../]
  [./topy]
    type = FunctionPresetBC
    boundary = top
    variable = disp_y
    function = pull
  [../]
  [./bottomt]
    type = PresetBC
    boundary = bottom
    variable = temp
    value = 300.0
  [../]
  [./topt]
    type = PresetBC
    boundary = top
    variable = temp
    value = 400.0
  [../]
[] # BCs

[Constraints]
  [./xfem_constraint]
    type = XFEMSingleVariableConstraint
    variable = temp
    jump = 0.0
    jump_flux = 0.0
    use_displaced_mesh = true
  [../]
[]


[Materials]
  [./elas]
    type = Elastic
    block = 0
    disp_x = disp_x
    disp_y = disp_y
    youngs_modulus = 75
    poissons_ratio = 0.3
    thermal_expansion = 26e-6
    formulation = PlaneStrain
    temp = temp
    stress_free_temperature = 300.0
  [../]
  [./heat]
    type = HeatConductionMaterial
    block = 0
    specific_heat = 1.0
    thermal_conductivity = 1.0
  [../]
[]

[Executioner]

  type = Transient

  solve_type = 'PJFNK'

  petsc_options = ksp_monitor
#  petsc_options_iname = '-pc_type -ksp_grmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
#  petsc_options_value = 'asm         31   preonly   lu      1'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -ksp_gmres_restart'
  petsc_options_value = 'lu       superlu_dist                  101'

  line_search = 'none'

   l_max_its = 50
   nl_max_its = 40

   nl_rel_step_tol= 1e-10
   nl_rel_tol = 1e-10


   start_time = 0.0
   dt = 1

   end_time = 1
   num_steps = 1

[]

[Outputs]
  exodus = true
[]

[Preconditioning]
  active = 'smp'
  [./smp]
    type = SMP
    pc_side = left
    full = true
  [../]
[]
