//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GrainDataTracker.h"
#include "RankTwoTensor.h"

class GrainTrackerRotationTensor;
class EulerAngleProvider;

template <>
InputParameters validParams<GrainTrackerRotationTensor>();

/**
 * Manage a list of elasticity tensors for the grains
 */
class GrainTrackerRotationTensor : public GrainDataTracker<RankTwoTensor>
{
public:
  GrainTrackerRotationTensor(const InputParameters & parameters);

protected:
  RankTwoTensor newGrain(unsigned int new_grain_id);

  /// generate random rotations when the Euler Angle provider runs out of data (otherwise error out)
  const bool _random_rotations;

  /// object providing the Euler angles
  const EulerAngleProvider & _euler;
};
