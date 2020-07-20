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
  InputParameters params = Postprocessor::validParams();
  //params += DiscreteElementUserObject::validParams();
  params.addClassDescription(
      "Retrieves the gradient value on the specified side of a specified interface ");
  params.addRequiredParam<UserObjectName>("value_at_interface_uo",
      "The name of the userobject that obtains the value and gradient at the interface.");
  params.addRequiredParam<Real>("side",
      "Side of the interface we want the gradient of (-1 for negative, +1 for positive)");
  return params;
}

GradValueAtXFEMInterfacePostprocessor::GradValueAtXFEMInterfacePostprocessor(const InputParameters & parameters)
  : //DiscreteElementUserObject(parameters),
    Postprocessor(parameters),
    _side(getParam<Real>("side"))
{
}

void
GradValueAtXFEMInterfacePostprocessor::initialize()
{
  _grad_value = 0;

  const UserObject * uo =
      &(_fe_problem.getUserObjectBase(getParam<UserObjectName>("value_at_interface_uo")));

  if (dynamic_cast<const PointValueAtXFEMInterface *>(uo) == nullptr)
    mooseError("UserObject casting to PointValueAtXFEMInterface in GradValueAtXFEMInterfacePostprocessor");

  _value_at_interface_uo = dynamic_cast<const PointValueAtXFEMInterface *>(uo);
}

void
GradValueAtXFEMInterfacePostprocessor::getGradientValue(unsigned int point_id)
{
  RealVectorValue grad_positive = _value_at_interface_uo->getGradientAtPositiveLevelSet()[point_id];
  RealVectorValue grad_negative = _value_at_interface_uo->getGradientAtNegativeLevelSet()[point_id];

  if (_side < 0)
  {
    _grad_value = grad_negative(0);
  }
  else
  {
    _grad_value = grad_positive(0);
  }
}

Real
GradValueAtXFEMInterfacePostprocessor::getValue()
{
  return _grad_value;
}
