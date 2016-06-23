/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef COMPUTEPOLYCRYSTALELASTICITYTENSORCP_H
#define COMPUTEPOLYCRYSTALELASTICITYTENSORCP_H

#include "ComputeElasticityTensorBase.h"

//Forward Declarations
class ComputePolycrystalElasticityTensorCP;
class GrainTracker;
class EulerAngleProvider;
class RotationTensor;

/**
 * Compute an evolving elasticity tensor coupled to a grain growth phase field model.
 */
class ComputePolycrystalElasticityTensorCP : public ComputeElasticityTensorBase
{
public:
  ComputePolycrystalElasticityTensorCP(const InputParameters & parameters);

protected:
  virtual void computeQpElasticityTensor();

  /// Unrotated elasticity tensor
  RankFourTensor _C_unrotated;

  /// Object providing the Euler angles
  const EulerAngleProvider & _euler;

  /// Grain tracker object
  const GrainTracker & _grain_tracker;

  /// Number of grains
  unsigned int _grain_num;

  /// Number of order parameters
  unsigned int _nop;

  /// Order parameters
  std::vector<const VariableValue *> _vals;

  MaterialProperty<std::vector<RankFourTensor> > & _elastic_tensor_cp;

  /// vector of elasticity tensors for storing rotated elasticity tensors for each grain
  std::vector<RankFourTensor> _C_rotated;

  std::vector<RankTwoTensor> _rot_trans;

  /// Crystal Rotation Matrix
  MaterialProperty<std::vector<RankTwoTensor> > & _crysrot;
};

#endif //COMPUTEPOLYCRYSTALELASTICITYTENSORCP_H
