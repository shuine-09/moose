//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GrainTrackerRotationTensor.h"
#include "EulerAngleProvider.h"
#include "RotationTensor.h"

registerMooseObject("PhaseFieldApp", GrainTrackerRotationTensor);

template <>
InputParameters
validParams<GrainTrackerRotationTensor>()
{
  InputParameters params = validParams<GrainTracker>();
  params.addParam<bool>("random_rotations",
                        true,
                        "Generate random rotations when the Euler Angle "
                        "provider runs out of data (otherwise error "
                        "out)");
  params.addRequiredParam<UserObjectName>("euler_angle_provider",
                                          "Name of Euler angle provider user object");
  return params;
}

GrainTrackerRotationTensor::GrainTrackerRotationTensor(const InputParameters & parameters)
  : GrainDataTracker<RankTwoTensor>(parameters),
    _random_rotations(getParam<bool>("random_rotations")),
    _euler(getUserObject<EulerAngleProvider>("euler_angle_provider"))
{
}

RankTwoTensor
GrainTrackerRotationTensor::newGrain(unsigned int new_grain_id)
{
  EulerAngles angles;

  if (new_grain_id < _euler.getGrainNum())
    angles = _euler.getEulerAngles(new_grain_id);
  else
  {
    if (_random_rotations)
      angles.random();
    else
      mooseError("GrainTrackerRotationTensor has run out of grain rotation data.");
  }

  return RotationTensor(RealVectorValue(angles));
}
