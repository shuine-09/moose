import h5py
dream3d = h5py.File('weld.dream3d', 'r')

euler_angles_file = open('euler_angles.txt', 'w+')
feature_ids_file = open('feature_ids.txt', 'w+')

#EulerAngles
euler_angles = dream3d['DataContainers']['SyntheticVolumeDataContainer']['CellData']['EulerAngles']

#FeatureIds
feature_ids = dream3d['DataContainers']['SyntheticVolumeDataContainer']['CellData']['FeatureIds']

#print(euler_angles.value.shape)

rg1 = (euler_angles.value.shape)[0]
rg2 = (euler_angles.value.shape)[1]
rg3 = (euler_angles.value.shape)[2]

for k in range(rg3):
  for j in range(rg2):
    for i in range(rg1):
      euler_angles_file.write("%f %f %f\n" % (euler_angles.value[i,j,k,0], euler_angles.value[i,j,k,1], euler_angles.value[i,j,k,2]))

dream3d_base = h5py.File('base.dream3d', 'r')

#EulerAngles
euler_angles_base = dream3d_base['DataContainers']['SyntheticVolumeDataContainer']['CellData']['EulerAngles']

#FeatureIds
feature_ids_base = dream3d_base['DataContainers']['SyntheticVolumeDataContainer']['CellData']['FeatureIds']

rg1 = (euler_angles_base.value.shape)[0]
rg2 = (euler_angles_base.value.shape)[1]
rg3 = (euler_angles_base.value.shape)[2]

for i in range(rg3):
  for j in range(rg2):
    for k in range(rg1, 2*rg1):
      euler_angles_file.write("%f %f %f\n" % (euler_angles_base.value[i,j,k-rg1,0], euler_angles_base.value[i,j,k-rg1,1], euler_angles_base.value[i,j,k-rg1,2]))

dream3d.close()
euler_angles_file.close()
feature_ids_file.close()
