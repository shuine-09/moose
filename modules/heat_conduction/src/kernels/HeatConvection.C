//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HeatConvection.h"

registerADMooseObject("HeatConductionApp", HeatConvection);

defineADLegacyParams(HeatConvection);

template <ComputeStage compute_stage>
InputParameters
HeatConvection<compute_stage>::validParams()
{
  InputParameters params = ADKernelGrad<compute_stage>::validParams();
  params.addRequiredCoupledVar("velocity", "The velocity");
  params.addParam<MaterialPropertyName>(
      "density", "rho", "the name of the density material property");
  params.addParam<MaterialPropertyName>(
      "enthalpy", "enthalpy", "the name of the enthalpy material property");
  params.addRequiredCoupledVar("level_set", "Level set variable");
  return params;
}

template <ComputeStage compute_stage>
HeatConvection<compute_stage>::HeatConvection(const InputParameters & parameters)
  : ADKernelGrad<compute_stage>(parameters),
    _velocity(adCoupledVectorValue("velocity")),
    _rho(getADMaterialProperty<Real>("density")),
    _h(getADMaterialProperty<Real>("enthalpy")),
    _ls(adCoupledValue("level_set"))
{
}

template <ComputeStage compute_stage>
ADRealVectorValue
HeatConvection<compute_stage>::precomputeQpResidual()
{
  // return -_rho[_qp] * _h[_qp] * _velocity[_qp] * std::abs(1 - _ls[_qp]);
  return -_rho[_qp] * _h[_qp] * _velocity[_qp];
}
