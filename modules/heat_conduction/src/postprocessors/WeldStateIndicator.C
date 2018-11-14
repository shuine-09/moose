//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "WeldStateIndicator.h"

#include <algorithm>
#include <limits>

registerMooseObject("MooseApp", WeldStateIndicator);

template <>
InputParameters
validParams<WeldStateIndicator>()
{
  // Define the parameters
  InputParameters params = validParams<ElementUserObject>();
  params.addRequiredCoupledVar("variable",
                               "The name of the variable that this postprocessor operates on");
  params.addRequiredParam<Real>("start_time", "Time to begin heating.");
  params.addRequiredParam<Real>("end_time", "Time to end heating.");
  params.addRequiredParam<Real>("melting_temp", "Melting temperature for weld material.");
  return params;
}

WeldStateIndicator::WeldStateIndicator(const InputParameters & parameters)
  : ElementUserObject(parameters),
    _u(coupledValue("variable")),
    _min_u(std::numeric_limits<Real>::max()),
    _has_melted(false),
    _start_time(getParam<Real>("start_time")),
    _end_time(getParam<Real>("end_time")),
    _melting_temp(getParam<Real>("melting_temp"))
{
}

void
WeldStateIndicator::initialize()
{
  _min_u = std::numeric_limits<Real>::max(); // start w/ the max
}

void
WeldStateIndicator::execute()
{
  for (unsigned int _qp = 0; _qp < _qrule->n_points(); _qp++)
    _min_u = std::min(_min_u, _u[_qp]);
}

void
WeldStateIndicator::threadJoin(const UserObject & y)
{
  const WeldStateIndicator & pps = static_cast<const WeldStateIndicator &>(y);

  _min_u = std::min(_min_u, pps._min_u);
}

void
WeldStateIndicator::finalize()
{
  gatherMin(_min_u);

  _min_u *= 1.5;

  if (_t < _start_time)
    _state = BEFORE;
  else if (_min_u < _melting_temp && _t < _end_time && !_has_melted)
    _state = HEATING;
  else if (_min_u >= _melting_temp && _t < _end_time)
  {
    _state = COOLING;
    _has_melted = true;
  }
  else
    _state = AFTER;
}
