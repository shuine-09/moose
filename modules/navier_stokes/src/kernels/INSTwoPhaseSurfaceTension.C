//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "INSTwoPhaseSurfaceTension.h"

registerADMooseObject("NavierStokesApp", INSTwoPhaseSurfaceTension);

defineADValidParams(INSTwoPhaseSurfaceTension,
                    ADVectorKernelValue,
                    params.addRequiredCoupledVar("grad_c",
                                                 "Regularized gradient of Level set variable");
                    params.addRequiredCoupledVar("curvature", "Level set variable");
                    params.addRequiredParam<Real>("sigma", "surface tension coefficient."););

template <ComputeStage compute_stage>
INSTwoPhaseSurfaceTension<compute_stage>::INSTwoPhaseSurfaceTension(
    const InputParameters & parameters)
  : ADVectorKernelValue<compute_stage>(parameters),
    _grad_c(adCoupledVectorValue("grad_c")),
    _curvature(adCoupledValue("curvature")),
    _sigma(getParam<Real>("sigma"))
{
}

template <ComputeStage compute_stage>
ADRealVectorValue
INSTwoPhaseSurfaceTension<compute_stage>::precomputeQpResidual()
{
  // ADReal curvature =
  //     -_second_c[_qp].tr() / _grad_c[_qp].norm() +
  //     _second_c[_qp] * _grad_c[_qp] * _grad_c[_qp] / std::pow(_grad_c[_qp].norm(), 3.0);
  //
  // return -_sigma * curvature * _grad_c[_qp];
  return -_sigma * _curvature[_qp] * _grad_c[_qp];
}
