/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "XFEMSingleVariableConstraintLS.h"
#include "FEProblem.h"
#include "DisplacedProblem.h"
#include "Assembly.h"
#include "XFEM.h"

// libMesh includes
#include "libmesh/quadrature.h"
#include "libmesh/fe_interface.h"

template<>
InputParameters validParams<XFEMSingleVariableConstraintLS>()
{
  InputParameters params = validParams<ElemElemConstraint>();
  params.addParam<Real>("alpha", 100, "Stablization parameter in Nitsche's formulation.");
  params.addParam<Real>("jump", 0, "Jump at the interface.");
  params.addParam<Real>("jump_flux", 0, "Flux jump at the interface.");
  params.addRequiredParam<NonlinearVariableName>("level_set_var", "The name of level set variable used to represent the interface");
  return params;
}

XFEMSingleVariableConstraintLS::XFEMSingleVariableConstraintLS(const InputParameters & parameters) :
    ElemElemConstraint(parameters),
    _alpha(getParam<Real>("alpha")),
    _jump(getParam<Real>("jump")),
    _jump_flux(getParam<Real>("jump_flux")),
    _ls_var_name(parameters.get<NonlinearVariableName>("level_set_var")),
    _ls_var(_fe_problem.getVariable(_tid, _ls_var_name)),
    _aux_system(_fe_problem.getAuxiliarySystem()),
    _aux_solution(_aux_system.currentSolution())
{
  FEProblem * fe_problem = dynamic_cast<FEProblem *>(&_subproblem);
  if (fe_problem == NULL)
    mooseError("Problem casting _subproblem to FEProblem in XFEMSingleVariableConstraintLS");
  _xfem = MooseSharedNamespace::dynamic_pointer_cast<XFEM>(fe_problem->getXFEM());
  if (_xfem == NULL)
    mooseError("Problem casting to XFEM in XFEMSingleVariableConstraintLS");
}

XFEMSingleVariableConstraintLS::~XFEMSingleVariableConstraintLS()
{
}

void
XFEMSingleVariableConstraintLS::reinitConstraintQuadrature(const ElementPairInfo & element_pair_info)
{
  _interface_normal = element_pair_info._elem1_normal;
  ElemElemConstraint::reinitConstraintQuadrature(element_pair_info);
}

bool
XFEMSingleVariableConstraintLS::shouldReverseSign()
{
  bool reverse = true;
  const Elem * undisplaced_elem  = NULL;
  FEProblem * _fe_problem = dynamic_cast<FEProblem *>(&_subproblem);
  if(_fe_problem->getDisplacedProblem() != NULL)
    undisplaced_elem = _fe_problem->getDisplacedProblem()->refMesh().elemPtr(_current_elem->id());
  else
    undisplaced_elem = _current_elem;

  for ( unsigned int i = 0; i < _current_elem->n_nodes(); ++i)
  {
    Real flag = _xfem->flagQpointInside(undisplaced_elem, *(_current_elem->get_node(i))); //inside (flag = 1) or ouside (flag = 0) real domain
    if (flag > 0.5)
    {
      dof_id_type ls_dof_id = (_current_elem->get_node(i))->dof_number(_aux_system.number(), _ls_var.number(), 0);
      Number ls_node = (*_aux_solution)(ls_dof_id);
      if (ls_node > 0.0)
        reverse = false;
      else
        reverse = true;

      break;
    }
  }
  return reverse;
}

Real
XFEMSingleVariableConstraintLS::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0;

  Real jump = _jump;
  Real jump_flux = _jump_flux;

  if  (shouldReverseSign())
  {
    jump *= -1.0;
    jump_flux *= -1.0;
  }

  switch (type)
  {
    case Moose::Element:
      r -= (0.5 * _grad_u[_qp] * _interface_normal + 0.5 * _grad_u_neighbor[_qp] * _interface_normal) * _test[_i][_qp];
      r -= (_u[_qp] - _u_neighbor[_qp]) * 0.5 * _grad_test[_i][_qp] * _interface_normal;
      r += 0.5 * _grad_test[_i][_qp] * _interface_normal * jump + 0.5 * _test[_i][_qp] * jump_flux;
      r += _alpha * (_u[_qp] - _u_neighbor[_qp] - jump) * _test[_i][_qp];
      break;

    case Moose::Neighbor:
      r += (0.5 * _grad_u[_qp] * _interface_normal + 0.5 * _grad_u_neighbor[_qp] * _interface_normal) * _test_neighbor[_i][_qp];
      r -= (_u[_qp] - _u_neighbor[_qp]) * 0.5 * _grad_test_neighbor[_i][_qp] * _interface_normal;
      r += 0.5 * _grad_test_neighbor[_i][_qp] * _interface_normal * jump + 0.5 * _test_neighbor[_i][_qp] * jump_flux;
      r -= _alpha * (_u[_qp] - _u_neighbor[_qp] - jump) * _test_neighbor[_i][_qp];
      break;
  }
  return r;
}

Real
XFEMSingleVariableConstraintLS::computeQpJacobian(Moose::DGJacobianType type)
{
  Real r = 0;

  switch (type)
  {
    case Moose::ElementElement:
      r += -0.5 * _grad_phi[_j][_qp] * _interface_normal * _test[_i][_qp] - _phi[_j][_qp] * 0.5 * _grad_test[_i][_qp] * _interface_normal;
      r += _alpha * _phi[_j][_qp] * _test[_i][_qp];
      break;

    case Moose::ElementNeighbor:
      r += -0.5 * _grad_phi_neighbor[_j][_qp] * _interface_normal * _test[_i][_qp] + _phi_neighbor[_j][_qp] * 0.5 * _grad_test[_i][_qp] * _interface_normal;
      r -= _alpha * _phi_neighbor[_j][_qp] * _test[_i][_qp];
      break;

    case Moose::NeighborElement:
      r += 0.5 * _grad_phi[_j][_qp] * _interface_normal * _test_neighbor[_i][_qp] - _phi[_j][_qp] * 0.5 * _grad_test_neighbor[_i][_qp] * _interface_normal;
      r -= _alpha * _phi[_j][_qp] * _test_neighbor[_i][_qp];
      break;

    case Moose::NeighborNeighbor:
      r += 0.5 * _grad_phi_neighbor[_j][_qp] * _interface_normal * _test_neighbor[_i][_qp] + _phi_neighbor[_j][_qp] * 0.5 * _grad_test_neighbor[_i][_qp] * _interface_normal;
      r += _alpha * _phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp];
      break;
  }

  return r;
}
