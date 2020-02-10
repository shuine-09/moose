//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// MOOSE includes
#include "LevelSetCurvatureRegulization.h"

registerADMooseObject("LevelSetApp", LevelSetCurvatureRegulization);

defineADValidParams(LevelSetCurvatureRegulization,
                    ADKernel,
                    params.addRequiredCoupledVar("grad_c", "The level set variable.");
                    params.addRequiredParam<Real>("epsilon", "The epsilon coefficient."););

template <ComputeStage compute_stage>
LevelSetCurvatureRegulization<compute_stage>::LevelSetCurvatureRegulization(
    const InputParameters & parameters)
  : ADKernel<compute_stage>(parameters),
    _grad_c(adCoupledVectorValue("grad_c")),
    _epsilon(getParam<Real>("epsilon"))
{
}

template <ComputeStage compute_stage>
ADReal
LevelSetCurvatureRegulization<compute_stage>::computeQpResidual()
{
  ADReal s = std::numeric_limits<ADReal>::epsilon();
  ADRealVectorValue n(0, 0, 0);
  if (_grad_c[_qp].norm_sq() > 1.0e-8)
  {
    s = _grad_c[_qp].norm();
    n = _grad_c[_qp] / s;
  }
  return _test[_i][_qp] * _u[_qp] - _grad_test[_i][_qp] * (n - _epsilon * _grad_u[_qp]);
}
