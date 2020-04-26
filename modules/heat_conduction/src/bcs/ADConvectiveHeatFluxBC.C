//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADConvectiveHeatFluxBC.h"

registerMooseObject("HeatConductionApp", ADConvectiveHeatFluxBC);

InputParameters
ADConvectiveHeatFluxBC::validParams()
{
  InputParameters params = ADIntegratedBC::validParams();
  params.addClassDescription(
      "Convective heat transfer boundary condition with temperature and heat "
      "transfer coefficent given by material properties.");
  params.addRequiredParam<Real>("heat_convection_coef", "Heat convection coeffient");
  params.addRequiredParam<Real>("ambient_temperature", "Ambient temperature");
  params.addRequiredParam<Real>("emissivity", "Emissivity");
  params.addRequiredParam<Real>("stefan_boltzmann", "Stefan boltzmann constant");
  params.addRequiredParam<Real>("substrate_thickness", "Substrate thickness");
  params.addRequiredParam<Real>("absorptivity", "Absorptivity");
  params.addRequiredParam<Real>("laser_power", "Laser power");
  params.addRequiredParam<Real>("laser_radius", "Laser radius");
  params.addRequiredParam<Real>("laser_velocity", "Laser velocity");
  return params;
}

ADConvectiveHeatFluxBC::ADConvectiveHeatFluxBC(const InputParameters & parameters)
  : ADIntegratedBC(parameters),
    _heat_convection_coef(getParam<Real>("heat_convection_coef")),
    _ambient_temperature(getParam<Real>("ambient_temperature")),
    _emissivity(getParam<Real>("emissivity")),
    _stefan_boltzmann(getParam<Real>("stefan_boltzmann")),
    _substrate_thickness(getParam<Real>("substrate_thickness")),
    _absorptivity(getParam<Real>("absorptivity")),
    _laser_power(getParam<Real>("laser_power")),
    _laser_radius(getParam<Real>("laser_radius")),
    _laser_velocity(getParam<Real>("laser_velocity"))
{
}

ADReal
ADConvectiveHeatFluxBC::computeQpResidual()
{
  ADReal I0 = 0;

  if (std::abs(_q_point[_qp](0) - _laser_velocity * _t) <= _laser_radius)
    I0 = _absorptivity * _laser_power / libMesh::pi / Utility ::pow<2>(_laser_radius);

  return -_test[_i][_qp] * (I0 - _heat_convection_coef * (_u[_qp] - _ambient_temperature) -
                            _emissivity * _stefan_boltzmann *
                                (Utility::pow<4>(_u[_qp]) - Utility::pow<4>(_ambient_temperature)));
}
