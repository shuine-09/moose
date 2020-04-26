//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADMatHeatSource.h"

registerMooseObject("HeatConductionApp", ADMatHeatSource);

InputParameters
ADMatHeatSource::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addRequiredParam<Real>("heat_convection_coef", "Heat convection coeffient");
  params.addRequiredParam<Real>("ambient_temperature", "Ambient temperature");
  params.addRequiredParam<Real>("emissivity", "Emissivity");
  params.addRequiredParam<Real>("stefan_boltzmann", "Stefan boltzmann constant");
  params.addRequiredParam<Real>("substrate_thickness", "Substrate thickness");
  return params;
}

ADMatHeatSource::ADMatHeatSource(const InputParameters & parameters)
  : ADKernel(parameters),
    _heat_convection_coef(getParam<Real>("heat_convection_coef")),
    _ambient_temperature(getParam<Real>("ambient_temperature")),
    _emissivity(getParam<Real>("emissivity")),
    _stefan_boltzmann(getParam<Real>("stefan_boltzmann")),
    _substrate_thickness(getParam<Real>("substrate_thickness"))
{
}

ADReal
ADMatHeatSource::computeQpResidual()
{
  return -_test[_i][_qp] *
         (-2.0 *
          (_heat_convection_coef * (_u[_qp] - _ambient_temperature) +
           _emissivity * _stefan_boltzmann *
               (Utility::pow<4>(_u[_qp]) - Utility::pow<4>(_ambient_temperature))) /
          _substrate_thickness);
}
