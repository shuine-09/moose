//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PFCrackVolume.h"

// MOOSE includes
#include "MooseVariable.h"

registerMooseObject("PhaseFieldApp", PFCrackVolume);

InputParameters
PFCrackVolume::validParams()
{
  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addRequiredParam<VariableName>("variable",
                                        "The name of the phase field fracture variable");
  params.addRequiredParam<Real>("l", "Internal length scale ");
  return params;
}

PFCrackVolume::PFCrackVolume(const InputParameters & parameters)
  : ElementIntegralPostprocessor(parameters),
    MooseVariableInterface<Real>(this,
                                 false,
                                 "variable",
                                 Moose::VarKindType::VAR_ANY,
                                 Moose::VarFieldType::VAR_FIELD_STANDARD),
    _var(_subproblem.getStandardVariable(_tid, parameters.get<VariableName>("variable"))),
    _u(_var.sln()),
    _grad_u(_var.gradSln()),
    _u_dot(_var.uDot()),
    _l(getParam<Real>("l")) // K
{
  addMooseVariableDependency(mooseVariable());
}

Real
PFCrackVolume::computeQpIntegral()
{
  return 1.0 / _l * _u[_qp] * _u[_qp] + _l * _grad_u[_qp] * _grad_u[_qp];
}
