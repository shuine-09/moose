//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DoubleEllipsoidHeatSource.h"

registerMooseObject("HeatConductionApp", DoubleEllipsoidHeatSource);

template <>
InputParameters
validParams<DoubleEllipsoidHeatSource>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<Real>("power", "laser power");
  params.addParam<Real>("efficienty", 1, "process efficienty");
  params.addRequiredParam<Real>("a", "transverse ellipsoid axe");
  params.addRequiredParam<Real>("b", "depth ellipsoid axe");
  params.addRequiredParam<Real>("c", "longitudinal ellipsoid axe");
  params.addRequiredParam<Real>("velocity", "heating spot travel speed");
  params.addParam<Real>("factor", 1, "scaling factor");

  params.addClassDescription("Double ellipsoid volumetric source heat.");

  return params;
}

DoubleEllipsoidHeatSource::DoubleEllipsoidHeatSource(const InputParameters & parameters)
  : Material(parameters),
    _P(getParam<Real>("power")),
    _eta(getParam<Real>("efficienty")),
    _a(getParam<Real>("a")),
    _b(getParam<Real>("b")),
    _c(getParam<Real>("c")),
    _v(getParam<Real>("velocity")),
    _f(getParam<Real>("factor")),
    _volumetric_heat(declareProperty<Real>("volumetric_heat"))
{
}

void
DoubleEllipsoidHeatSource::computeQpProperties()
{
  const Real & x = _q_point[_qp](0);
  const Real & y = _q_point[_qp](1);
  const Real & z = _q_point[_qp](2);
  _volumetric_heat[_qp] = 6.0 * std::sqrt(3.0) * _P * _eta * _f /
                          (_a * _b * _c * std::pow(libMesh::pi, 1.5)) *
                          std::exp(-(3.0 * std::pow(x, 2.0) / std::pow(_a, 2.0) +
                                     3.0 * std::pow(y, 2.0) / std::pow(_b, 2.0) +
                                     3.0 * std::pow(z + _v * _t, 2.0) / pow(_c, 2.0)));
}
