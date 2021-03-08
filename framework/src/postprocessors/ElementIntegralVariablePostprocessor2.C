//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElementIntegralVariablePostprocessor2.h"

registerMooseObject("MooseApp", ElementIntegralVariablePostprocessor2);

defineLegacyParams(ElementIntegralVariablePostprocessor2);

InputParameters
ElementIntegralVariablePostprocessor2::validParams()
{
  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addRequiredCoupledVar("variable", "The name of the variable that this object operates on");
  params.addClassDescription("Computes a volume integral of the specified variable");
  params.addRequiredCoupledVar("displacements",
                               "The string of displacements suitable for the problem statement");
  return params;
}

ElementIntegralVariablePostprocessor2::ElementIntegralVariablePostprocessor2(
    const InputParameters & parameters)
  : ElementIntegralPostprocessor(parameters),
    MooseVariableInterface<Real>(this,
                                 false,
                                 "variable",
                                 Moose::VarKindType::VAR_ANY,
                                 Moose::VarFieldType::VAR_FIELD_STANDARD),
    _u(coupledValue("variable")),
    _grad_u(coupledGradient("variable")),
    _ndisp(coupledComponents("displacements")),
    _grad_disp(3),
    _disp(3),
    _disp_var(3),
    _pressure(getMaterialProperty<Real>("fracture_pressure"))
{
  addMooseVariableDependency(mooseVariable());
  for (unsigned int i = 0; i < _ndisp; ++i)
  {
    _grad_disp[i] = &coupledGradient("displacements", i);
    _disp[i] = &coupledValue("displacements", i);
    _disp_var[i] = coupled("displacements", i);
  }
  for (unsigned i = _ndisp; i < 3; ++i)
  {
    _disp[i] = &_zero;
    _grad_disp[i] = &_grad_zero;
  }
}

Real
ElementIntegralVariablePostprocessor2::computeQpIntegral()
{
  return _pressure[_qp] * _grad_u[_qp].norm() * 2.0 * _u[_qp];
}
