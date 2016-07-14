/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "XFEMDisplacementConstraint.h"
#include "FEProblem.h"
#include "Assembly.h"

// libMesh includes
#include "libmesh/quadrature.h"

template<>
InputParameters validParams<XFEMDisplacementConstraint>()
{
  InputParameters params = validParams<ElemElemConstraint>();
  params.addParam<Real>("alpha", 100, "Stablization parameter in Nitsche's formulation.");
  return params;
}

XFEMDisplacementConstraint::XFEMDisplacementConstraint(const InputParameters & parameters) :
    ElemElemConstraint(parameters),
    _alpha(getParam<Real>("alpha"))
{
}

XFEMDisplacementConstraint::~XFEMDisplacementConstraint()
{
}

void
XFEMDisplacementConstraint::reinitConstraintQuadrature(const ElementPairInfo & element_pair_info)
{
  _interface_normal = element_pair_info._elem1_normal;
  ElemElemConstraint::reinitConstraintQuadrature(element_pair_info);
}

Real
XFEMDisplacementConstraint::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0;
 
  //std::cout << "u size = " << _u.size() << std::endl;
  //std::cout << "u_neighbor size = " << _u_neighbor.size() << std::endl;

  //std::cout << "u[" << _qp << "] " << _u[_qp] << std::endl;
  //std::cout << "u_neighbor[ " << _qp  << "] = " << _u_neighbor[_qp] << std::endl;

  switch (type)
  {
    case Moose::Element:
      r += _alpha * (_u[_qp] - _u_neighbor[_qp]) * _test[_i][_qp];
      break;

    case Moose::Neighbor:
      r -= _alpha * (_u[_qp] - _u_neighbor[_qp]) * _test_neighbor[_i][_qp];
      break;
  }

  //std::cout << "r = " << r << std::endl;
  return r;
}

Real
XFEMDisplacementConstraint::computeQpJacobian(Moose::DGJacobianType type)
{
  Real r = 0;

  switch (type)
  {
    case Moose::ElementElement:
      r += _alpha * _phi[_j][_qp] * _test[_i][_qp];
      break;

    case Moose::ElementNeighbor:
      r -= _alpha * _phi_neighbor[_j][_qp] * _test[_i][_qp];
      break;

    case Moose::NeighborElement:
      r -= _alpha * _phi[_j][_qp] * _test_neighbor[_i][_qp];
      break;

    case Moose::NeighborNeighbor:
      r += _alpha * _phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp];
      break;
  }

  return r;
}
