//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// MOOSE includes
#include "LevelSetPowderMass.h"

registerADMooseObject("LevelSetApp", LevelSetPowderMass);

defineADValidParams(LevelSetPowderMass,
                    ADKernelValue,
                    params.addClassDescription("Implement powder material mass addition.");
                    params.addParam<Real>("mass_rate", 2.5e-4, "Mass rate.");
                    params.addParam<Real>("mass_radius", 0.25e-3, "Mass radius.");
                    params.addRequiredParam<UserObjectName>("laser_location",
                                                            "Userobject name of laser location."););

template <ComputeStage compute_stage>
LevelSetPowderMass<compute_stage>::LevelSetPowderMass(const InputParameters & parameters)
  : ADKernelValue<compute_stage>(parameters),
    _rho(getADMaterialProperty<Real>("rho")),
    _mass_rate(getParam<Real>("mass_rate")),
    _mass_radius(getParam<Real>("mass_radius")),
    _location(
        getUserObjectByName<MeltPoolLevelSetLocation>(getParam<UserObjectName>("laser_location")))
{
}

template <ComputeStage compute_stage>
ADReal
LevelSetPowderMass<compute_stage>::precomputeQpResidual()
{
  ADRealVectorValue laser_center = _location.getLaserSpotLocation();
  ADReal r = (_ad_q_point[_qp] - laser_center).norm();
  ADReal power_feed = 0;

  if (r <= _mass_radius)
    power_feed = _mass_rate * std::exp(-Utility::pow<2>(r / _mass_radius));

  return (_grad_u[_qp] + RealVectorValue(1.0e-10)).norm() * power_feed;
}
