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
    _volumetric_heat(declareADProperty<Real>("volumetric_heat"))
{
}

void
DoubleEllipsoidHeatSource::computeQpProperties()
{
  // const Real & x = _q_point[_qp](0);
  // const Real & y = _q_point[_qp](1);
  // // const Real & z = _q_point[_qp](2);
  // Real x_t = 22;
  // Real y_t = 10;
  // if (_v * _t > 0 && _v * _t < 50)
  // {
  //   x_t = 22 + _v * _t;
  //   y_t = 10;
  // }
  // else if (_v * _t > 50 && _v * _t < 100)
  // {
  //   x_t = 78 - (_v * _t - 50);
  //   y_t = 12;
  // }
  // else if (_v * _t > 100 && _v * _t < 150)
  // {
  //   x_t = 22 + (_v * _t - 100);
  //   y_t = 14;
  // }
  // else if (_v * _t > 150 && _v * _t < 200)
  // {
  //   x_t = 78 - (_v * _t - 150);
  //   y_t = 16;
  // }
  // else if (_v * _t > 200 && _v * _t < 250)
  // {
  //   x_t = 22 + (_v * _t - 200);
  //   y_t = 18;
  // }
  // else if (_v * _t > 250 && _v * _t < 300)
  // {
  //   x_t = 78 - (_v * _t - 250);
  //   y_t = 20;
  // }
  // else if (_v * _t > 300 && _v * _t < 350)
  // {
  //   x_t = 22 + (_v * _t - 300);
  //   y_t = 22;
  // }
  // else if (_v * _t > 350 && _v * _t < 400)
  // {
  //   x_t = 78 - (_v * _t - 350);
  //   y_t = 24;
  // }
  // else if (_v * _t > 400 && _v * _t < 450)
  // {
  //   x_t = 22 + (_v * _t - 400);
  //   y_t = 26;
  // }
  // else if (_v * _t > 450 && _v * _t < 500)
  // {
  //   x_t = 78 - (_v * _t - 450);
  //   y_t = 28;
  // }
  // else if (_v * _t > 500 && _v * _t < 550)
  // {
  //   x_t = 22 + (_v * _t - 500);
  //   y_t = 30;
  // }
  // else
  // {
  //   x_t = 1000;
  //   y_t = 1000;
  // }
  //
  // _volumetric_heat[_qp] = 6.0 * std::sqrt(3.0) * _P * _eta * _f /
  //                         (_a * _b * _c * std::pow(libMesh::pi, 1.5)) *
  //                         std::exp(-(3.0 * std::pow(x - x_t, 2.0) / std::pow(_a, 2.0) +
  //                                    3.0 * std::pow(y - y_t, 2.0) / std::pow(_b, 2.0) +
  //                                    3.0 * std::pow(0, 2.0) / pow(_c, 2.0)));
  const Real & x = _q_point[_qp](0);
  const Real & y = _q_point[_qp](1);
  const Real & z = _q_point[_qp](2);

  Real x_t = 9 * std::cos(2 * libMesh::pi * _t);
  Real y_t = 9 * std::sin(2 * libMesh::pi * _t);
  Real z_t = std::floor(_t) * 1 + 5;

  _volumetric_heat[_qp] = 6.0 * std::sqrt(3.0) * _P * _eta * _f /
                          (_a * _b * _c * std::pow(libMesh::pi, 1.5)) *
                          std::exp(-(3.0 * std::pow(x - x_t, 2.0) / std::pow(_a, 2.0) +
                                     3.0 * std::pow(y - y_t, 2.0) / std::pow(_b, 2.0) +
                                     3.0 * std::pow(z - z_t, 2.0) / std::pow(_c, 2.0)));
}
