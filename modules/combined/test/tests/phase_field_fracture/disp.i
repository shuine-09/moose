[Mesh]
  file = fuel_3d.e
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  block = '1'
[]

[Variables]
  [./temp]
  [../]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./disp_z]
  [../]
[]

[AuxVariables]
  [stress_radial]
    family = MONOMIAL
    order = CONSTANT
  []
  [stress_hoop]
    family = MONOMIAL
    order = CONSTANT
  []
  [stress_axial]
    family = MONOMIAL
    order = CONSTANT
  []
  [./gc]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./d]
  [../]
  [./hist]
   order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
   [./gc]
     type = MaterialRealAux
     variable = gc
     property = gc_prop
   [../]
  [./radial]
    type = RankTwoScalarAux
    variable = stress_radial
    rank_two_tensor = stress
    scalar_type = RadialStress
    point1 = '0 0 0'
    point2 = '0 1 0'
    execute_on = timestep_end
  [../]
  [./hoop]
    type = RankTwoScalarAux
    variable = stress_hoop
    rank_two_tensor = stress
    scalar_type = HoopStress
    point1 = '0 0 0'
    point2 = '0 1 0'
    execute_on = timestep_end
  [../]
  [./axial]
    type = RankTwoScalarAux
    variable = stress_axial
    rank_two_tensor = stress
    scalar_type = AxialStress
    point1 = '0 0 0'
    point2 = '0 1 0'
    execute_on = timestep_end
  [../]
  [./hist]
    type = MaterialRealAux
    variable = hist
    property = hist
    execute_on = timestep_end
  [../]

[]

[Kernels]
  [Hcond]
    type = HeatConduction
    variable = temp
    coupled_to_damage = false
    c = d
  []
  [Hsource]
    type = HeatSource
    variable = temp
    function = power_history
  []
[]


[Modules]
  [./TensorMechanics]
    [./Master]
      [./mech]
        add_variables = false
        strain = SMALL
        additional_generate_output = 'stress_yy stress_xx stress_zz strain_xx strain_yy strain_zz vonmises_stress'
        volumetric_locking_correction = false
        #decomposition_method = EigenSolution
        eigenstrain_names = thermal_eigenstrain
        temperature = temp
      [../]
    [../]
  [../]
[]

[BCs]
  [no_x]
    type = DirichletBC
    variable = disp_x
    boundary = '100'
    value = 0.0
  []
  [no_y]
    type = DirichletBC
    variable = disp_y
    boundary = '100'
    value = 0.0
  []
  [no_z]
    type = DirichletBC
    variable = disp_z
    boundary = '100'
    value = 0.0
  []

#  [no_x]
#    type = DirichletBC
#    variable = disp_x
#    boundary = '4 61'
#    value = 0.0
#  []
#  [no_z]
#    type = DirichletBC
#    variable = disp_z
#    boundary = '4 51'
#    value = 0.0
#  []
#  [no_y]
#    type = DirichletBC
#    variable = disp_y
#    boundary = '4'
#    value = 0.0
#  []
  [pellet_outer_temp]
    type = FunctionDirichletBC
    variable = temp
    boundary = 1
    function = temp_or
  []
[]

[MultiApps]
  [./sub]
    type = TransientMultiApp
    execute_on = timestep_end
    positions = '0 0 0' ## how how to define this
    input_files = pf.i
  [../]
[]

[Transfers]
  [./from_sub]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    multi_app = sub
    source_variable = d
    variable = d
  [../]
  [./to_sub_x]
    type = MultiAppCopyTransfer
    direction = to_multiapp
    multi_app = sub
    source_variable = disp_x
    variable = disp_x
  [../]
  [./to_sub_y]
    type = MultiAppCopyTransfer
    direction = to_multiapp
    multi_app = sub
    source_variable = disp_y
    variable = disp_y
  [../]
  [./to_sub_z]
    type = MultiAppCopyTransfer
    direction = to_multiapp
    multi_app = sub
    source_variable = disp_z
    variable = disp_z
  [../]
  [./to_sub_t]
    type = MultiAppCopyTransfer
    direction = to_multiapp
    multi_app = sub
    source_variable = temp
    variable = temp
  [../]
[]

[Materials]
#   [./gc]
#    type = WeibullMaterial
#    weibull_modulus = 15
#    specimen_volume = 1.25e-13
#    specimen_material_property = 40  #18.2
#  [../]
  [./pfbulkmat]
    type = GenericConstantMaterial
    prop_names = 'l visco'
    prop_values = '3e-4 100'
  [../]
  [gc]
    type = GenericFunctionMaterial
    prop_values = 'gc'
    prop_names = 'gc_prop'
  []

   [./elastic]
    type = ComputeIsotropicLinearElasticPFFractureStress
    c = d
    F_name = F
    use_current_history_variable = false
    kdamage = 1e-8
    decomp = true
    linear = true
    use_vi = true
   [../]
#   [./elastic]
#     type = ComputeLinearElasticPFFractureStress
#     c = d
#     F_name = F
#     use_current_history_variable = false
#     decomp = false
#     kdamage = 1e-5
#  [../]
  [./define_mobility]
    type = ParsedMaterial
    material_property_names = 'visco gc_prop'
    f_name = L
    function = '1.0/(visco * gc_prop)'
  [../]
  [./define_kappa]
    type = ParsedMaterial
    material_property_names = 'gc_prop l'
    f_name = kappa_op
    function = 'gc_prop * l'
  [../]

  [fuel_density]
    type = Density
    block = '1'
    density = 10431.0
  []
  [fuel_thermal]
    type = HeatConductionMaterial
    block = '1'
    temp = 'temp'
    thermal_conductivity = 5.0
    specific_heat = 1.0
  []
  [Cijkl]
    type = ComputeIsotropicElasticityTensor
    poissons_ratio = 0.3
    youngs_modulus = 2e11
  []
  [thermal_expansion]
    type = ComputeThermalExpansionEigenstrain
    stress_free_temperature = '300'
    temperature = 'temp'
    thermal_expansion_coeff = 10e-6
    eigenstrain_name = thermal_eigenstrain
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    ksp_norm = default
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  # petsc_options_iname = '-ksp_type -pc_type -pc_factor_mat_solver_package'
  # petsc_options_value = 'preonly   lu       superlu_dist'
 #petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
#petsc_options_value = 'lu       superlu_dist'

  petsc_options_iname = '-pc_type -sub_pc_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold'
  petsc_options_value = 'hypre    boomeramg    201                0.9'

  picard_max_its = 10
  picard_rel_tol = 1e-8
  picard_abs_tol = 1e-8

  nl_abs_tol = 1e-8
  nl_rel_tol = 1e-8
  l_max_its = 30
  nl_max_its = 8
  dt = 5
  dtmax = 5
  end_time = 10000
  line_search = none
  # [./TimeStepper]
  #   type = IterationAdaptiveDT
  #   dt = 5
  #   optimal_iterations = 5
  #   iteration_window = 0
  #   growth_factor = 1.5
  #   cutback_factor = 0.5
  # [../]
[]

[Outputs]
  print_linear_residuals = false
  [./exodus]
    type = Exodus
    interval = 1
  [../]
  [./my_checkpoint]
    type = Checkpoint
    num_files = 2
    interval = 100
  [../]
  csv = true
[]

[Functions]
 [gc0]
    type = PiecewiseMultilinear
    data_file = 'gc.txt'
  []
  [gc]
    type = CompositeFunction
    functions = 'gc0'
    scale_factor = 1
  []
  [power_history]
    type = PiecewiseLinear
    scale_factor = 1.8936e+04 # Convert from W/m to W/m3 (1/pi/(.0041m)^2/0.0115)
    x = '0   1e4    1e8'
    y = '0   2.5e4  2.5e4'
  []
  [temp_or]
    type = PiecewiseLinear
    xy_data = '0.0000000E+000   3.0000000E+002
               1.0000000E+002   5.7602620E+002
               2.0000000E+002   5.8868136E+002
               3.0000000E+002   5.9354501E+002
               4.0000000E+002   5.9803698E+002
               5.0000000E+002   6.0246166E+002
               6.0000000E+002   6.0683317E+002
               7.0000000E+002   6.1115311E+002
               8.0000000E+002   6.1542251E+002
               9.0000000E+002   6.1964236E+002
               1.0000000E+003   6.2381358E+002
               1.1000000E+003   6.2793709E+002
               1.2000000E+003   6.3201376E+002
               1.3000000E+003   6.3604445E+002
               1.4000000E+003   6.4002998E+002
               1.5000000E+003   6.4397114E+002
               1.6000000E+003   6.4786872E+002
               1.7000000E+003   6.5172346E+002
               1.8000000E+003   6.5553611E+002
               1.9000000E+003   6.5930737E+002
               2.0000000E+003   6.6303796E+002
               2.1000000E+003   6.6672853E+002
               2.2000000E+003   6.7037976E+002
               2.3000000E+003   6.7399234E+002
               2.4000000E+003   6.7756681E+002
               2.5000000E+003   6.8110388E+002
               2.6000000E+003   6.8460411E+002
               2.7000000E+003   6.8806812E+002
               2.8000000E+003   6.9149648E+002
               2.9000000E+003   6.9488976E+002
               3.0000000E+003   6.9824854E+002
               3.1000000E+003   7.0157338E+002
               3.2000000E+003   7.0486476E+002
               3.3000000E+003   7.0812322E+002
               3.4000000E+003   7.1134929E+002
               3.5000000E+003   7.1454345E+002
               3.6000000E+003   7.1770619E+002
               3.7000000E+003   7.2083798E+002
               3.8000000E+003   7.2393919E+002
               3.9000000E+003   7.2701042E+002
               4.0000000E+003   7.3005203E+002
               4.1000000E+003   7.3306446E+002
               4.2000000E+003   7.3604808E+002
               4.3000000E+003   7.3900325E+002
               4.4000000E+003   7.4193037E+002
               4.5000000E+003   7.4482980E+002
               4.6000000E+003   7.4770188E+002
               4.7000000E+003   7.5054697E+002
               4.8000000E+003   7.5336536E+002
               4.9000000E+003   7.5615738E+002
               5.0000000E+003   7.5892332E+002
               5.1000000E+003   7.6166351E+002
               5.2000000E+003   7.6437809E+002
               5.3000000E+003   7.6706738E+002
               5.4000000E+003   7.6973168E+002
               5.5000000E+003   7.7237121E+002
               5.6000000E+003   7.7498620E+002
               5.7000000E+003   7.7757687E+002
               5.8000000E+003   7.8014344E+002
               5.9000000E+003   7.8268613E+002
               6.0000000E+003   7.8520510E+002
               6.1000000E+003   7.8770060E+002
               6.2000000E+003   7.9017276E+002
               6.3000000E+003   7.9262179E+002
               6.4000000E+003   7.9504785E+002
               6.5000000E+003   7.9745112E+002
               6.6000000E+003   7.9983176E+002
               6.7000000E+003   8.0218986E+002
               6.8000000E+003   8.0452571E+002
               6.9000000E+003   8.0683937E+002
               7.0000000E+003   8.0913095E+002
               7.1000000E+003   8.1140067E+002
               7.2000000E+003   8.1364861E+002
               7.3000000E+003   8.1587485E+002
               7.4000000E+003   8.1807958E+002
               7.5000000E+003   8.2026287E+002
               7.6000000E+003   8.2242480E+002
               7.7000000E+003   8.2456544E+002
               7.8000000E+003   8.2668483E+002
               7.9000000E+003   8.2047428E+002
               8.0000000E+003   8.1649166E+002
               8.1000000E+003   8.1519224E+002
               8.2000000E+003   8.1470729E+002
               8.3000000E+003   8.1459804E+002
               8.4000000E+003   8.1471836E+002
               8.5000000E+003   8.1499084E+002
               8.6000000E+003   8.1536815E+002
               8.7000000E+003   8.1582063E+002
               8.8000000E+003   8.1632844E+002
               8.9000000E+003   8.1687609E+002
               9.0000000E+003   8.1745282E+002
               9.1000000E+003   8.1805006E+002
               9.2000000E+003   8.1866094E+002
               9.3000000E+003   8.1927988E+002
               9.4000000E+003   8.1990224E+002
               9.5000000E+003   8.2052362E+002
               9.6000000E+003   8.2114187E+002
               9.7000000E+003   8.2175388E+002
               9.8000000E+003   8.2235750E+002
               9.9000000E+003   8.2295103E+002
               1.0000000E+004   8.2353311E+002
               1.0100000E+004   8.2350435E+002
               1.0200000E+004   8.2344182E+002
               1.0300000E+004   8.2338339E+002
               1.0400000E+004   8.2332938E+002
               1.0500000E+004   8.2327904E+002
               1.0600000E+004   8.2323176E+002
               1.0700000E+004   8.2318705E+002
               1.0800000E+004   8.2314457E+002
               1.0900000E+004   8.2310404E+002
               1.1000000E+004   8.2306521E+002
               1.1100000E+004   8.2302790E+002
               1.1200000E+004   8.2299193E+002
               1.1300000E+004   8.2295719E+002
               1.1400000E+004   8.2292354E+002
               1.1500000E+004   8.2289090E+002
               1.1600000E+004   8.2285919E+002
               1.1700000E+004   8.2282831E+002
               1.1800000E+004   8.2279822E+002
               1.1900000E+004   8.2276885E+002
               1.2000000E+004   8.2274015E+002
               1.2100000E+004   8.2271209E+002
               1.2200000E+004   8.2268460E+002
               1.2300000E+004   8.2265767E+002
               1.2400000E+004   8.2263126E+002
               1.2500000E+004   8.2260534E+002
               1.2600000E+004   8.2257988E+002
               1.2700000E+004   8.2255485E+002
               1.2800000E+004   8.2253024E+002
               1.2900000E+004   8.2250602E+002
               1.3000000E+004   8.2248217E+002
               1.3100000E+004   8.2245868E+002
               1.3200000E+004   8.2243553E+002
               1.3300000E+004   8.2241270E+002
               1.3400000E+004   8.2239019E+002
               1.3500000E+004   8.2236797E+002
               1.3600000E+004   8.2234604E+002
               1.3700000E+004   8.2232439E+002
               1.3800000E+004   8.2230300E+002
               1.3900000E+004   8.2228186E+002
               1.4000000E+004   8.2226096E+002
               1.4100000E+004   8.2224031E+002
               1.4200000E+004   8.2221987E+002
               1.4300000E+004   8.2219966E+002
               1.4400000E+004   8.2217966E+002
               1.4500000E+004   8.2215986E+002
               1.4600000E+004   8.2214027E+002
               1.4700000E+004   8.2212086E+002
               1.4800000E+004   8.2210165E+002
               1.4900000E+004   8.2208262E+002
               1.5000000E+004   8.2206375E+002'
  []
[]
