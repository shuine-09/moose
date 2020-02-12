//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MeltPoolHeatSource.h"
#include "Function.h"

registerADMooseObject("HeatConductionApp", MeltPoolHeatSource);

defineADLegacyParams(MeltPoolHeatSource);

template <ComputeStage compute_stage>
InputParameters
MeltPoolHeatSource<compute_stage>::validParams()
{
  InputParameters params = ADKernel<compute_stage>::validParams();
  params.addRequiredCoupledVar("level_set", "Level set variable");
  params.addRequiredParam<UserObjectName>("laser_location", "Userobject name of laser location.");
  params.addRequiredParam<Real>("laser_power", "Laser power.");
  params.addRequiredParam<Real>("effective_beam_radius", "Effective beam radius.");
  params.addRequiredParam<Real>("absorption_coefficient", "Absorption coefficient.");
  params.addRequiredParam<Real>("heat_transfer_coefficient", "Heat transfer coefficient.");
  params.addRequiredParam<Real>("StefanBoltzmann_constant", "Stefan Boltzmann constant.");
  params.addRequiredParam<Real>("material_emissivity", "Material emissivity.");
  params.addRequiredParam<Real>("ambient_temperature", "Ambient temperature.");
  return params;
}

template <ComputeStage compute_stage>
MeltPoolHeatSource<compute_stage>::MeltPoolHeatSource(const InputParameters & parameters)
  : ADKernel<compute_stage>(parameters),
    _grad_ls(adCoupledGradient("level_set")),
    _power(getParam<Real>("laser_power")),
    _alpha(getParam<Real>("absorption_coefficient")),
    _Rb(getParam<Real>("effective_beam_radius")),
    _Ah(getParam<Real>("heat_transfer_coefficient")),
    _stefan_boltzmann(getParam<Real>("StefanBoltzmann_constant")),
    _varepsilon(getParam<Real>("material_emissivity")),
    _T0(getParam<Real>("ambient_temperature")),
    _location(
        getUserObjectByName<MeltPoolLevelSetLocation>(getParam<UserObjectName>("laser_location")))
{
}

template <ComputeStage compute_stage>
ADReal
MeltPoolHeatSource<compute_stage>::computeQpResidual()
{
  ADRealVectorValue laser_center = _location.getLaserSpotLocation();
  ADReal r = (_ad_q_point[_qp] - laser_center).norm();

  ADReal laser_source = 2 * _power * _alpha / (libMesh::pi * Utility::pow<2>(_Rb)) *
                        std::exp(-2.0 * Utility::pow<2>(r / _Rb));

  ADReal convection = -_Ah * (_u[_qp] - _T0);
  ADReal radiation =
      -_stefan_boltzmann * _varepsilon * (Utility::pow<4>(_u[_qp]) - Utility::pow<4>(_T0));

  ADReal heat_source =
      (convection + radiation + laser_source) * (_grad_ls[_qp] + RealVectorValue(1.0e-10)).norm();
  return -_test[_i][_qp] * heat_source;
}
