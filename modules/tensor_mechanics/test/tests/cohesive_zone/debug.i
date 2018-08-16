[Mesh]
  file = n10-id1.msh
# file = grains_debug.e
#  file = poly_3d.e
#  file = gmsh_test_in.e
#  file = hex_grains_v2.e
  parallel_type = REPLICATED
[]

[MeshModifiers]
  [./breakmesh]
    type = BreakMeshByBlock
    split_interface = true
  [../]
[]

