//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VariableScaleAux.h"

registerMooseObject("MooseApp", VariableScaleAux);

template <>
InputParameters
validParams<VariableScaleAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("Scale varibale based on a level set function");
  params.addRequiredParam<VariableName>("scale_var", "The variable from which to scale");
  params.addParam<Real>("scale_factor", 1, "Scale factor.");
  return params;
}

VariableScaleAux::VariableScaleAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _scale_var_number(_subproblem
                          .getVariable(_tid,
                                       parameters.get<VariableName>("scale_var"),
                                       Moose::VarKindType::VAR_ANY,
                                       Moose::VarFieldType::VAR_FIELD_STANDARD)
                          .number()),
    _system(_subproblem.getSystem(getParam<VariableName>("scale_var"))),
    _solution(_system.current_local_solution.get()),
    _scale_factor(getParam<Real>("scale_factor"))
{
  if (!isNodal())
    mooseError("VariableScaleAux variable is NOT nodal");
}

Real
VariableScaleAux::computeValue()
{
  const Node * node = _current_node;

  dof_id_type dof_id = node->dof_number(_system.number(), _scale_var_number, 0);
  Number node_value = (*_solution)(dof_id);

  return node_value * _scale_factor;
}
