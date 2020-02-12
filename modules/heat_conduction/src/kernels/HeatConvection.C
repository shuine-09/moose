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
      "specific_heat", "specific_heat", "the name of the specific heat material property");
  return params;
}

template <ComputeStage compute_stage>
HeatConvection<compute_stage>::HeatConvection(const InputParameters & parameters)
  : ADKernelGrad<compute_stage>(parameters),
    _velocity(adCoupledVectorValue("velocity")),
    _rho(getADMaterialProperty<Real>("density")),
    _cp(getADMaterialProperty<Real>("specific_heat"))
{
}

template <ComputeStage compute_stage>
ADRealVectorValue
HeatConvection<compute_stage>::precomputeQpResidual()
{
  return -_rho[_qp] * _cp[_qp] * _velocity[_qp];
}
