//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// MOOSE includes
#include "LevelSetOlssonPlane.h"

registerMooseObject("LevelSetApp", LevelSetOlssonPlane);

template <>
InputParameters
validParams<LevelSetOlssonPlane>()
{
  InputParameters params = validParams<Function>();
  params.addClassDescription("Implementation of 'bubble' ranging from 0 to 1.");
  params.addParam<Real>("epsilon", 0.01, "The interface thickness.");
  params.addParam<bool>("z", false, "The interface thickness.");
  return params;
}

LevelSetOlssonPlane::LevelSetOlssonPlane(const InputParameters & parameters)
  : Function(parameters), _epsilon(getParam<Real>("epsilon")), _z(getParam<bool>("z"))
{
}

Real
LevelSetOlssonPlane::value(Real /*t*/, const Point & p) const
{
  Real x = -(p(1) - 0.005) / _epsilon;
  if (_z)
    x = -(p(2) - 0.005) / _epsilon;
  return 1.0 / (1 + std::exp(x));
}

RealGradient
LevelSetOlssonPlane::gradient(Real /*t*/, const Point & p) const
{
  return p;
}
