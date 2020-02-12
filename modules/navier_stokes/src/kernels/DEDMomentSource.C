//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DEDMomentSource.h"

registerADMooseObject("NavierStokesApp", DEDMomentSource);

defineADValidParams(DEDMomentSource,
                    ADVectorKernelValue,
                    params.addRequiredParam<MaterialPropertyName>("ded_momentum", "DED momentum"););

template <ComputeStage compute_stage>
DEDMomentSource<compute_stage>::DEDMomentSource(const InputParameters & parameters)
  : ADVectorKernelValue<compute_stage>(parameters),
    _ded_momentum(getADMaterialProperty<RealVectorValue>("ded_momentum"))
{
}

template <ComputeStage compute_stage>
ADRealVectorValue
DEDMomentSource<compute_stage>::precomputeQpResidual()
{
  return -_ded_momentum[_qp];
}
