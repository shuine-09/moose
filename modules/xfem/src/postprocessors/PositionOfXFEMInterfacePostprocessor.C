//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PositionOfXFEMInterfacePostprocessor.h"

registerMooseObject("XFEMApp", PositionOfXFEMInterfacePostprocessor);

InputParameters
PositionOfXFEMInterfacePostprocessor::validParams()
{
  InputParameters params = GeneralPostprocessor::validParams();
  params.addClassDescription(
      "Retrieves the gradient value on the specified side of a specified interface ");
  params.addRequiredParam<UserObjectName>("value_at_interface_uo",
      "The name of the userobject that obtains the value and gradient at the interface.");
  return params;
}

PositionOfXFEMInterfacePostprocessor::PositionOfXFEMInterfacePostprocessor(const InputParameters & parameters)
  : GeneralPostprocessor(parameters)
{
}

void
PositionOfXFEMInterfacePostprocessor::initialize()
{
  const UserObject * uo =
      &(_fe_problem.getUserObjectBase(getParam<UserObjectName>("value_at_interface_uo")));

  if (dynamic_cast<const PointValueAtXFEMInterface *>(uo) == nullptr)
    mooseError("UserObject casting to PointValueAtXFEMInterface in PositionOfXFEMInterfacePostprocessor");

  _value_at_interface_uo = dynamic_cast<const PointValueAtXFEMInterface *>(uo);
}


Real
PositionOfXFEMInterfacePostprocessor::getValue()
{
  return _value_at_interface_uo->getCurrentX();
}
