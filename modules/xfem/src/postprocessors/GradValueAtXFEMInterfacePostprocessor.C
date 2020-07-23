//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GradValueAtXFEMInterfacePostprocessor.h"

registerMooseObject("XFEMApp", GradValueAtXFEMInterfacePostprocessor);

InputParameters
GradValueAtXFEMInterfacePostprocessor::validParams()
{
  InputParameters params = GeneralPostprocessor::validParams();
  params.addClassDescription(
      "Retrieves the gradient value on the specified side of a specified interface ");
  params.addRequiredParam<UserObjectName>("value_at_interface_uo",
      "The name of the userobject that obtains the value and gradient at the interface.");
  params.addRequiredParam<Real>("side",
      "Side of the interface we want the gradient of (-1 for negative, +1 for positive)");
  return params;
}

GradValueAtXFEMInterfacePostprocessor::GradValueAtXFEMInterfacePostprocessor(const InputParameters & parameters)
  : GeneralPostprocessor(parameters),
    _side(getParam<Real>("side"))
{
}

void
GradValueAtXFEMInterfacePostprocessor::initialize()
{
  const UserObject * uo =
      &(_fe_problem.getUserObjectBase(getParam<UserObjectName>("value_at_interface_uo")));

  if (dynamic_cast<const PointValueAtXFEMInterface *>(uo) == nullptr)
    mooseError("UserObject casting to PointValueAtXFEMInterface in GradValueAtXFEMInterfacePostprocessor");

  _value_at_interface_uo = dynamic_cast<const PointValueAtXFEMInterface *>(uo);
}


Real
GradValueAtXFEMInterfacePostprocessor::getValue()
{
  if (_side < 0)
  {
    //std::cout << "Negative side" << std::endl;
    return _value_at_interface_uo->getGradientXComponentAtNegativeLevelSet();
  }
  else
  {
    //std::cout << "Positive side" << std::endl;
    return _value_at_interface_uo->getGradientXComponentAtPositiveLevelSet();
  }
}
