/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "XFEMSingleVariableConstraintPts.h"
#include "FEProblem.h"
#include "Assembly.h"
#include "XFEM.h"

// libMesh includes
#include "libmesh/quadrature.h"

template<>
InputParameters validParams<XFEMSingleVariableConstraintPts>()
{
  InputParameters params = validParams<ElemElemConstraint>();
  params.addParam<Real>("alpha", 100, "Stablization parameter in Nitsche's formulation.");
  params.addParam<Real>("jump", 0, "Jump at the interface.");
  params.addParam<Real>("jump_flux", 0, "Flux jump at the interface.");
  return params;
}

XFEMSingleVariableConstraintPts::XFEMSingleVariableConstraintPts(const InputParameters & parameters) :
    ElemElemConstraint(parameters),
    _alpha(getParam<Real>("alpha")),
    _jump(getParam<Real>("jump")),
    _jump_flux(getParam<Real>("jump_flux"))
{
  FEProblem * fe_problem = dynamic_cast<FEProblem *>(&_subproblem);
  if (fe_problem == NULL)
    mooseError("Problem casting _subproblem to FEProblem in XFEMSingleVariableConstraintPts");
  _xfem = MooseSharedNamespace::dynamic_pointer_cast<XFEM>(fe_problem->getXFEM());
  if (_xfem == NULL)
    mooseError("Problem casting to XFEM in XFEMSingleVariableConstraintPts");

}

XFEMSingleVariableConstraintPts::~XFEMSingleVariableConstraintPts()
{
}

void
XFEMSingleVariableConstraintPts::reinitConstraintQuadrature(const ElementPairInfo & element_pair_info)
{
  _interface_normal = element_pair_info._elem1_normal;
  ElemElemConstraint::reinitConstraintQuadrature(element_pair_info);
}

Real
XFEMSingleVariableConstraintPts::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0;

  switch (type)
  {
    case Moose::Element:
      r -= (0.5 * _grad_u[_qp] * _interface_normal + 0.5 * _grad_u_neighbor[_qp] * _interface_normal) * _test[_i][_qp];
      r -= (_u[_qp] - _u_neighbor[_qp]) * 0.5 * _grad_test[_i][_qp] * _interface_normal;
      r += 0.5 * _grad_test[_i][_qp] * _interface_normal * _jump + 0.5 * _test[_i][_qp] * _jump_flux;
      r += _alpha * (_u[_qp] - _u_neighbor[_qp] - _jump) * _test[_i][_qp];
      break;

    case Moose::Neighbor:
      r += (0.5 * _grad_u[_qp] * _interface_normal + 0.5 * _grad_u_neighbor[_qp] * _interface_normal) * _test_neighbor[_i][_qp];
      r -= (_u[_qp] - _u_neighbor[_qp]) * 0.5 * _grad_test_neighbor[_i][_qp] * _interface_normal;
      r += 0.5 * _grad_test_neighbor[_i][_qp] * _interface_normal * _jump + 0.5 * _test_neighbor[_i][_qp] * _jump_flux;
      r -= _alpha * (_u[_qp] - _u_neighbor[_qp] - _jump) * _test_neighbor[_i][_qp];
      break;
  }

  for(unsigned int i = 0; i < _xfem->numberInterfacePoints(); ++i)
  {
    Point p = _xfem->getInterfacePoint(i);
    if (_current_elem->contains_point(p))
    {
      _xfem->setInterfaceQuantity(i, (_u[0] + _u[1]) * 0.5); 
    }
  }
  return r;
}

Real
XFEMSingleVariableConstraintPts::computeQpJacobian(Moose::DGJacobianType type)
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
