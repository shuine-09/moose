//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LevelSetGradientRegulization.h"

registerADMooseObject("LevelSetApp", LevelSetGradientRegulization);

defineADLegacyParams(LevelSetGradientRegulization);

template <ComputeStage compute_stage>
InputParameters
LevelSetGradientRegulization<compute_stage>::validParams()
{
  InputParameters params = ADVectorKernel<compute_stage>::validParams();
  params.addRequiredCoupledVar("c", "Level set variable");
  params.addParam<bool>("normal_regulization", false, "Level set variable");
  return params;
}

template <ComputeStage compute_stage>
LevelSetGradientRegulization<compute_stage>::LevelSetGradientRegulization(
    const InputParameters & parameters)
  : ADVectorKernel<compute_stage>(parameters),
    _grad_c(adCoupledGradient("c")),
    _normal_regulization(getParam<bool>("normal_regulization"))
{
}

template <ComputeStage compute_stage>
ADReal
LevelSetGradientRegulization<compute_stage>::computeQpResidual()
{
  if (!_normal_regulization)
    return _test[_i][_qp] * (_u[_qp] - _grad_c[_qp]);
  else
  {
    if (_grad_c[_qp].norm_sq() > 1.0e-8)
      return _test[_i][_qp] * (_u[_qp] - _grad_c[_qp] / _grad_c[_qp].norm());
    else
      return _test[_i][_qp] * _u[_qp];
  }
}

template class LevelSetGradientRegulization<RESIDUAL>;
template class LevelSetGradientRegulization<JACOBIAN>;
