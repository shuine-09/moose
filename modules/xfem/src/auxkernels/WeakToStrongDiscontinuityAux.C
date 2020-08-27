//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "WeakToStrongDiscontinuityC4ZrAux.h"
#include "AuxiliarySystem.h"
#include "MooseVariable.h"
#include "XFEM.h"
#include "FEProblem.h"

registerMooseObject("MooseApp", WeakToStrongDiscontinuityC4ZrAux);

/**
 * Computes the value of the original variable of the strong discontinuity
 * problem from the weak discontinuity variable.
 * !!! Not working correctly. Only adds the discontinuity value once an element
 * has been fully crossed by the interface. !!!

 * Currently adds the value everywhere ! (to be used with a clip filter in
 * combination with athe weak variable. Only works for 1 inerface pb.)
 */

InputParameters
WeakToStrongDiscontinuityC4ZrAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription(
      "Auxiliary Kernel that calcuates the original problem variable by adding the discontinuity value on one side of the interface");
  params.addRequiredParam<VariableName>(
      "weak_variable", "The name of the variable solved for in the weak discontinuity problem");
  params.addRequiredParam<VariableName>(
      "level_set_var", "The name of level set variable used to represent the interface");
  params.addRequiredParam<Real>("jump_value","The value of the discontinuity");
  return params;
}

WeakToStrongDiscontinuityC4ZrAux::WeakToStrongDiscontinuityC4ZrAux(const InputParameters & parameters)
  : AuxKernel(parameters),
  _jump_value(getParam<Real>("jump_value")),
  _weak_variable_number(
      _subproblem.getVariable(_tid, parameters.get<VariableName>("weak_variable")).number()),
  _level_set_var_number(
      _subproblem.getVariable(_tid, parameters.get<VariableName>("level_set_var")).number()),
  _system(_subproblem.getSystem(getParam<VariableName>("level_set_var"))),
  _system2(_subproblem.getSystem(getParam<VariableName>("weak_variable"))),
  _solution(_system.current_local_solution.get()),
  _solution2(_system2.current_local_solution.get()),
  _add_jump(false)
{
  FEProblemBase * fe_problem = dynamic_cast<FEProblemBase *>(&_subproblem);

  _xfem = MooseSharedNamespace::dynamic_pointer_cast<XFEM>(fe_problem->getXFEM());
}

Real
WeakToStrongDiscontinuityC4ZrAux::computeValue()
{
  const Node * node = _current_node;

  dof_id_type ls_dof_id = node->dof_number(_system.number(), _level_set_var_number, 0);
  Number ls_node_value = (*_solution)(ls_dof_id);

  dof_id_type wv_dof_id = node->dof_number(_system2.number(), _weak_variable_number, 0);
  Number wv_node_value = (*_solution2)(wv_dof_id);

  //std::cout << "wv_node_value = " << wv_node_value << std::endl;
  //std::cout << "ls_node_value = " << ls_node_value << std::endl;

  _add_jump = false;

  if (ls_node_value > 0.0)
    //std::cout << "Physical node in the oxide" << std::endl;
    _add_jump = true;

  //std::cout << "_add_jump =" << _add_jump << std::endl;

  //return wv_node_value + _jump_value;


  if (_add_jump)
  {
    //std::cout << "u_discontinuous value in the oxide :" << wv_node_value + _jump_value << std::endl;
    return wv_node_value + _jump_value;
  }
  else
    return wv_node_value;
}
