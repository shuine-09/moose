//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MeltPoolMomentumSource.h"

registerADMooseObject("NavierStokesApp", MeltPoolMomentumSource);

defineADValidParams(MeltPoolMomentumSource,
                    ADVectorKernelValue,
                    params.addClassDescription("Melt pool momentum source term kernel."););

template <ComputeStage compute_stage>
MeltPoolMomentumSource<compute_stage>::MeltPoolMomentumSource(const InputParameters & parameters)
  : ADVectorKernelValue<compute_stage>(parameters),
    _melt_pool_momentum_source(getADMaterialProperty<RealVectorValue>("melt_pool_momentum_source"))
{
}

template <ComputeStage compute_stage>
ADRealVectorValue
MeltPoolMomentumSource<compute_stage>::precomputeQpResidual()
{
  return -_melt_pool_momentum_source[_qp];
}
