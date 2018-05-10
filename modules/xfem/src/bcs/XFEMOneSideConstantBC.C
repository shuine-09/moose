//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "XFEMOneSideConstantBC.h"

registerMooseObject("XFEMApp", XFEMOneSideConstantBC);

template <>
InputParameters
validParams<XFEMOneSideConstantBC>()
{
  InputParameters p = validParams<PresetBC>();
  p.addRequiredParam<Real>("left_x", "The left boundary");
  p.set<bool>("use_displaced_mesh") = false;
  return p;
}

XFEMOneSideConstantBC::XFEMOneSideConstantBC(const InputParameters & parameters)
  : PresetBC(parameters), _left_x(getParam<Real>("left_x"))
{
}

bool
XFEMOneSideConstantBC::shouldApply()
{
  if ((*_current_node)(0) < _left_x && _t_step == 1)
    return true;
  else
    return false;
}
