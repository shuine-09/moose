//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "RandomGc.h"

registerMooseObject("TensorMechanicsApp", RandomGc);

template <>
InputParameters
validParams<RandomGc>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<Real>("gc_mean", "mean gc property");
  params.addParam<Real>(
      "random_range", 0.0, "Range of a uniform random distribution for the threshold");
  return params;
}

RandomGc::RandomGc(const InputParameters & parameters)
  : Material(parameters),
    _gc(declareProperty<Real>("gc_prop")),
    _gc_old(getMaterialPropertyOld<Real>("gc_prop")),
    _gc_mean(getParam<Real>("gc_mean")),
    _random_range(getParam<Real>("random_range"))
{
  setRandomResetFrequency(EXEC_INITIAL);
}

void
RandomGc::initQpStatefulProperties()
{
}

void
RandomGc::computeQpProperties()
{
  if (_t_step == 1)
    _gc[_qp] = _gc_mean * ((1.0 - _random_range / 2.0) + _random_range * getRandomReal());
  else
    _gc[_qp] = _gc_old[_qp];
}
