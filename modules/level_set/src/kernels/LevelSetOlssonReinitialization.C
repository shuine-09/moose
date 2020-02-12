//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// MOOSE includes
#include "LevelSetOlssonReinitialization.h"

registerADMooseObject("LevelSetApp", LevelSetOlssonReinitialization);

defineADValidParams(
    LevelSetOlssonReinitialization,
    ADKernelGrad,
    params.addClassDescription("The re-initialization equation defined by Olsson et. al. (2007).");
    params.addRequiredCoupledVar(
        "phi_0", "The level set variable to be reinitialized as signed distance function.");
    params.addRequiredCoupledVar("vel_x", "vel_x.");
    params.addRequiredCoupledVar("vel_y", "vel_y.");
    params.addRequiredParam<PostprocessorName>("time",
                                               "The postprocessor associated with the C1111 value");
    params.addRequiredParam<PostprocessorName>(
        "epsilon", "The epsilon coefficient to be used in the reinitialization calculation.");
    params.addRequiredCoupledVar("grad_c", "The level set variable."););

template <ComputeStage compute_stage>
LevelSetOlssonReinitialization<compute_stage>::LevelSetOlssonReinitialization(
    const InputParameters & parameters)
  : ADKernelGrad<compute_stage>(parameters),
    _grad_levelset_0(adCoupledGradient("phi_0")),
    _levelset_0(adCoupledValue("phi_0")),
    _epsilon(getPostprocessorValue("epsilon")),
    _grad_c(adCoupledVectorValue("grad_c")),
    _vel_x(coupledValue("vel_x")),
    _vel_y(coupledValue("vel_y")),
    _master_time(getPostprocessorValue("time"))
{
}

template <ComputeStage compute_stage>
ADRealVectorValue
LevelSetOlssonReinitialization<compute_stage>::precomputeQpResidual()
{
  unsigned int step = std::floor(_master_time / 1.0e-4);

  // ADReal s = _grad_levelset_0[_qp].norm() + std::numeric_limits<ADReal>::epsilon();
  // ADRealVectorValue n_hat = _grad_levelset_0[_qp] / s;
  ADRealVectorValue n_hat = _grad_c[_qp];
  ADRealVectorValue f = _u[_qp] * (1 - _u[_qp]) * n_hat;
  // Real vel_mag = std::sqrt(_vel_x[_qp] * _vel_x[_qp] + _vel_y[_qp] * _vel_y[_qp]);
  // if (_levelset_0[_qp] > 0.99 || _levelset_0[_qp] < 0.01)
  //   return 0;
  // else if (vel_mag < 1.0e-8 && ((step % 10) != 0))
  //   return 0;
  // else
  return -f + _epsilon * (_grad_u[_qp]);
  // return (-f + _epsilon * (_grad_u[_qp] * n_hat) * n_hat);
}
