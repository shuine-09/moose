[GlobalParams]
  disp_x = disp_x
  disp_y = disp_y
  order = FIRST
  family = LAGRANGE
[]

[Mesh]
  file = unstructure_6.e
  displacements = 'disp_x disp_y'
  #uniform_refine = 1
[]

[XFEM]
  qrule = volfrac
  output_cut_plane = true
  use_crack_growth_increment = true
  crack_growth_increment = 0.0005 # old : 0.0005
[]

[UserObjects]
  [./xfem_marker_uo]
    type = XFEMMaterialTensorMarkerUserObject
    execute_on = timestep_end
    tensor = stress
    #quantity = MaxPrincipal
    quantity = hoop
    point1 = '0, 0, 0'
    point2 = '0, 0, 1'
    #average = false
    average = true
    threshold = 1.5e8 #1.3e8
    initiate_on_boundary = 8
    random_range = 0.1
    #weibull = xfem_weibull
    use_weibull = true
  [../]
  [./xfem_max_hoop_stress]
    type = XFEMMaxHoopStress
    execute_on = timestep_end
    block = pellet_type_1
    disp_x = disp_x
    disp_y = disp_y
    temp = temp
    youngs_modulus = 2.e11
    poissons_ratio = .345
    thermal_expansion = 10e-6
    radius_inner = 2e-4
    radius_outer = 4e-4
    average_h = h
    intersecting_boundary = 8
    critical_k = 8e6
 [../]
[]

[Variables]
  [./disp_x]
  [../]

  [./disp_y]
  [../]

  [./temp]
    initial_condition = 300.0     # set initial temp to ambient
  [../]
[]

[AuxVariables]
  [./xfem_volfrac]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./xfem_marker]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./xfem_cut_origin_x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./xfem_cut_origin_y]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./xfem_cut_origin_z]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./xfem_cut_normal_x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./xfem_cut_normal_y]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./xfem_cut_normal_z]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./xfem_threshold]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xx]      # stress aux variables are defined for output; this is a way to get integration point variables to the output file
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./hoop]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./radial]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./axial]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vonmises]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./weibull]
    order = CONSTANT
    family = MONOMIAL
  [../]

[]

[Functions]
  [./power_history]
    type = PiecewiseLinear
    scale_factor = 1.8936e+4 #Convert from W/m to W/m3 (1/pi/(.0041m)^2)
    x = '0,  1e4,    1e8'
    y = '0, 2.5e4, 2.5e4'
  [../]

  [./temp_or]
    type = PiecewiseLinear
    xy_data = '0.0000000E+000,  3.0000000E+002
               1.0000000E+002,  5.7602620E+002
               2.0000000E+002,  5.8868136E+002
               3.0000000E+002,  5.9354501E+002
               4.0000000E+002,  5.9803698E+002
               5.0000000E+002,  6.0246166E+002
               6.0000000E+002,  6.0683317E+002
               7.0000000E+002,  6.1115311E+002
               8.0000000E+002,  6.1542251E+002
               9.0000000E+002,  6.1964236E+002
               1.0000000E+003,  6.2381358E+002
               1.1000000E+003,  6.2793709E+002
               1.2000000E+003,  6.3201376E+002
               1.3000000E+003,  6.3604445E+002
               1.4000000E+003,  6.4002998E+002
               1.5000000E+003,  6.4397114E+002
               1.6000000E+003,  6.4786872E+002
               1.7000000E+003,  6.5172346E+002
               1.8000000E+003,  6.5553611E+002
               1.9000000E+003,  6.5930737E+002
               2.0000000E+003,  6.6303796E+002
               2.1000000E+003,  6.6672853E+002
               2.2000000E+003,  6.7037976E+002
               2.3000000E+003,  6.7399234E+002
               2.4000000E+003,  6.7756681E+002
               2.5000000E+003,  6.8110388E+002
               2.6000000E+003,  6.8460411E+002
               2.7000000E+003,  6.8806812E+002
               2.8000000E+003,  6.9149648E+002
               2.9000000E+003,  6.9488976E+002
               3.0000000E+003,  6.9824854E+002
               3.1000000E+003,  7.0157338E+002
               3.2000000E+003,  7.0486476E+002
               3.3000000E+003,  7.0812322E+002
               3.4000000E+003,  7.1134929E+002
               3.5000000E+003,  7.1454345E+002
               3.6000000E+003,  7.1770619E+002
               3.7000000E+003,  7.2083798E+002
               3.8000000E+003,  7.2393919E+002
               3.9000000E+003,  7.2701042E+002
               4.0000000E+003,  7.3005203E+002
               4.1000000E+003,  7.3306446E+002
               4.2000000E+003,  7.3604808E+002
               4.3000000E+003,  7.3900325E+002
               4.4000000E+003,  7.4193037E+002
               4.5000000E+003,  7.4482980E+002
               4.6000000E+003,  7.4770188E+002
               4.7000000E+003,  7.5054697E+002
               4.8000000E+003,  7.5336536E+002
               4.9000000E+003,  7.5615738E+002
               5.0000000E+003,  7.5892332E+002
               5.1000000E+003,  7.6166351E+002
               5.2000000E+003,  7.6437809E+002
               5.3000000E+003,  7.6706738E+002
               5.4000000E+003,  7.6973168E+002
               5.5000000E+003,  7.7237121E+002
               5.6000000E+003,  7.7498620E+002
               5.7000000E+003,  7.7757687E+002
               5.8000000E+003,  7.8014344E+002
               5.9000000E+003,  7.8268613E+002
               6.0000000E+003,  7.8520510E+002
               6.1000000E+003,  7.8770060E+002
               6.2000000E+003,  7.9017276E+002
               6.3000000E+003,  7.9262179E+002
               6.4000000E+003,  7.9504785E+002
               6.5000000E+003,  7.9745112E+002
               6.6000000E+003,  7.9983176E+002
               6.7000000E+003,  8.0218986E+002
               6.8000000E+003,  8.0452571E+002
               6.9000000E+003,  8.0683937E+002
               7.0000000E+003,  8.0913095E+002
               7.1000000E+003,  8.1140067E+002
               7.2000000E+003,  8.1364861E+002
               7.3000000E+003,  8.1587485E+002
               7.4000000E+003,  8.1807958E+002
               7.5000000E+003,  8.2026287E+002
               7.6000000E+003,  8.2242480E+002
               7.7000000E+003,  8.2456544E+002
               7.8000000E+003,  8.2668483E+002
               7.9000000E+003,  8.2047428E+002
               8.0000000E+003,  8.1649166E+002
               8.1000000E+003,  8.1519224E+002
               8.2000000E+003,  8.1470729E+002
               8.3000000E+003,  8.1459804E+002
               8.4000000E+003,  8.1471836E+002
               8.5000000E+003,  8.1499084E+002
               8.6000000E+003,  8.1536815E+002
               8.7000000E+003,  8.1582063E+002
               8.8000000E+003,  8.1632844E+002
               8.9000000E+003,  8.1687609E+002
               9.0000000E+003,  8.1745282E+002
               9.1000000E+003,  8.1805006E+002
               9.2000000E+003,  8.1866094E+002
               9.3000000E+003,  8.1927988E+002
               9.4000000E+003,  8.1990224E+002
               9.5000000E+003,  8.2052362E+002
               9.6000000E+003,  8.2114187E+002
               9.7000000E+003,  8.2175388E+002
               9.8000000E+003,  8.2235750E+002
               9.9000000E+003,  8.2295103E+002
               1.0000000E+004,  8.2353311E+002
               1.0100000E+004,  8.2350435E+002
               1.0200000E+004,  8.2344182E+002
               1.0300000E+004,  8.2338339E+002
               1.0400000E+004,  8.2332938E+002
               1.0500000E+004,  8.2327904E+002
               1.0600000E+004,  8.2323176E+002
               1.0700000E+004,  8.2318705E+002
               1.0800000E+004,  8.2314457E+002
               1.0900000E+004,  8.2310404E+002
               1.1000000E+004,  8.2306521E+002
               1.1100000E+004,  8.2302790E+002
               1.1200000E+004,  8.2299193E+002
               1.1300000E+004,  8.2295719E+002
               1.1400000E+004,  8.2292354E+002
               1.1500000E+004,  8.2289090E+002
               1.1600000E+004,  8.2285919E+002
               1.1700000E+004,  8.2282831E+002
               1.1800000E+004,  8.2279822E+002
               1.1900000E+004,  8.2276885E+002
               1.2000000E+004,  8.2274015E+002
               1.2100000E+004,  8.2271209E+002
               1.2200000E+004,  8.2268460E+002
               1.2300000E+004,  8.2265767E+002
               1.2400000E+004,  8.2263126E+002
               1.2500000E+004,  8.2260534E+002
               1.2600000E+004,  8.2257988E+002
               1.2700000E+004,  8.2255485E+002
               1.2800000E+004,  8.2253024E+002
               1.2900000E+004,  8.2250602E+002
               1.3000000E+004,  8.2248217E+002
               1.3100000E+004,  8.2245868E+002
               1.3200000E+004,  8.2243553E+002
               1.3300000E+004,  8.2241270E+002
               1.3400000E+004,  8.2239019E+002
               1.3500000E+004,  8.2236797E+002
               1.3600000E+004,  8.2234604E+002
               1.3700000E+004,  8.2232439E+002
               1.3800000E+004,  8.2230300E+002
               1.3900000E+004,  8.2228186E+002
               1.4000000E+004,  8.2226096E+002
               1.4100000E+004,  8.2224031E+002
               1.4200000E+004,  8.2221987E+002
               1.4300000E+004,  8.2219966E+002
               1.4400000E+004,  8.2217966E+002
               1.4500000E+004,  8.2215986E+002
               1.4600000E+004,  8.2214027E+002
               1.4700000E+004,  8.2212086E+002
               1.4800000E+004,  8.2210165E+002
               1.4900000E+004,  8.2208262E+002
               1.5000000E+004,  8.2206375E+002'
  [../]
[]

[SolidMechanics]
  [./solid]
    temp = temp
    volumetric_locking_correction = false
  [../]
[]

[Kernels]
  [./heat]         # gradient term in heat conduction equation
    type = HeatConduction
    variable = temp
  [../]

#  [./heat_ie]       # time term in heat conduction equation
#    type = HeatConductionTimeDerivative
#    variable = temp
#  [../]

  [./heat_source]  # source term in heat conduction equation
     type = HeatSource
     variable = temp
     block = pellet_type_1
     function = power_history
  [../]
[]

[AuxKernels]
 [./weibull_output]
   type = MaterialRealAux
   variable = weibull
   property = weibull
 [../]
  [./xfem_marker]
    type = XFEMMarkerAux
    variable = xfem_marker
    execute_on = linear
  [../]
 [./xfem_threshold]
    type = MaterialTensorAux
    variable = xfem_threshold
    tensor = stress
    quantity = MaxPrincipal
    #quantity = hoop
    #point1 = '0, 0, 0'
    #point2 = '0, 0, 1'
    execute_on = timestep_end
  [../]
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
  [./stress_zz]
    type = MaterialTensorAux
    tensor = stress
    variable = stress_zz
    index = 2
    execute_on = timestep_end
  [../]
  [./hoop]
    type = MaterialTensorAux
    tensor = stress
    variable = hoop
    quantity = hoop
    point1 = '0, 0, 0'
    point2 = '0, 0, 1'
    execute_on = linear
  [../]
  [./radial]
    type = MaterialTensorAux
    tensor = stress
    variable = radial
    quantity = radial
    point1 = '0, 0, 0'
    point2 = '0, 0, 1'
    execute_on = timestep_end
  [../]
  [./axial]
    type = MaterialTensorAux
    tensor = stress
    variable = axial
    quantity = axial
    point1 = '0, 0, 0'
    point2 = '0, 0, 1'
    execute_on = timestep_end
  [../]
  [./vonmises]
    type = MaterialTensorAux
    tensor = stress
    variable = vonmises
    quantity = vonmises
    execute_on = timestep_end
  [../]
[]

[BCs]
# Define boundary conditions

  [./no_x_center]
    type = DirichletBC
    variable = disp_x
    boundary = 17
    value = 0.0
  [../]

  [./no_y_center]
    type = DirichletBC
    variable = disp_y
    boundary = 17
    value = 0.0
  [../]

  [./no_y]
    type = DirichletBC
    variable = disp_y
    boundary = 18
    value = 0.0
  [../]

  [./pellet_outer_temp]
    type = FunctionPresetBC
    variable = temp
    boundary = 8
    function = temp_or
  [../]
[]

[Materials]
  [./fuel_thermal]
    type = HeatConductionMaterial
    block = pellet_type_1
    temp = temp
    thermal_conductivity = 5.0
    specific_heat = 1.0
  [../]

  [./fuel_elastic]
    type = Elastic
    block = pellet_type_1
    temp = temp
    youngs_modulus = 2.e11
    poissons_ratio = .345
    thermal_expansion = 10e-6
    compute_JIntegral = true
    stress_free_temperature = 300
  [../]

  [./fuel_density]
    type = Density
    block = pellet_type_1
    density = 10431.0
  [../]

  [./weilbull]
    type = XFEMWeibullMaterial
    weibull_modulus = 15
    specimen_volume = 5.281e-5 #5.281e-5 # 7.6977e-05
    specimen_material_property = 1.0
    block = pellet_type_1
  [../]
[]

[Dampers]
  [./limitT]
    type = MaxIncrement
    max_increment = 100.0
    variable = temp
  [../]
  [./limitX]
    type = MaxIncrement
    max_increment = 1e-5
    variable = disp_x
  [../]
[]

[Executioner]

  type = Transient

  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu superlu_dist'

  line_search = 'none'

  l_max_its = 100
  l_tol = 8e-3

  nl_max_its = 30
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-10

  start_time = 0.0
  dt = 1.0e2
  end_time = 1.5e4
  num_steps = 5000

  dtmin = 1.0e2

  max_xfem_update = 1
[]

[Postprocessors]
  [./_dt]                     # time step
    type = TimestepSize
  [../]

  [./nonlinear_its]           # number of nonlinear iterations at each timestep
    type = NumNonlinearIterations
  [../]

  [./h]
    type = AverageElementSize
    variable = temp
  [../]
[]

[Outputs]
  # Define output file(s)
  file_base = crack_init_5mm_grow_1mm
  interval = 1
  exodus = true
  [./console]
    type = Console
    perf_log = true
  [../]
[]
