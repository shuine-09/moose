//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// MOOSE includes
#include "LevelSetDEDMass.h"

registerADMooseObject("LevelSetApp", LevelSetDEDMass);

defineADValidParams(LevelSetDEDMass,
                    ADKernelValue,
                    params.addClassDescription("Implement powder material mass addition."););

template <ComputeStage compute_stage>
LevelSetDEDMass<compute_stage>::LevelSetDEDMass(const InputParameters & parameters)
  : ADKernelValue<compute_stage>(parameters),
    _mass_change(getADMaterialProperty<Real>("mass_change")),
    _powder_feed(getADMaterialProperty<Real>("powder_feed")),
    _rho(getADMaterialProperty<Real>("rho"))
{
}

template <ComputeStage compute_stage>
ADReal
LevelSetDEDMass<compute_stage>::precomputeQpResidual()
{
  return (_grad_u[_qp] + RealVectorValue(1.0e-10)).norm() * _powder_feed[_qp];
}
