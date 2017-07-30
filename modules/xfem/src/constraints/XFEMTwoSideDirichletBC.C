/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "XFEMTwoSideDirichletBC.h"

// MOOSE includes
#include "Assembly.h"
#include "ElementPairInfo.h"
#include "FEProblem.h"

// libMesh includes
#include "libmesh/quadrature.h"

template <>
InputParameters
validParams<XFEMTwoSideDirichletBC>()
{
  InputParameters params = validParams<ElemElemConstraint>();
  params.addParam<Real>("alpha", 100, "Penalty coefficient");
  params.addParam<Real>("value1", 0, "Dirichlet BC for one side");
  params.addParam<Real>("value2", 0, "Dirichlet BC for the other side");
  return params;
}

XFEMTwoSideDirichletBC::XFEMTwoSideDirichletBC(const InputParameters & parameters)
  : ElemElemConstraint(parameters),
    _alpha(getParam<Real>("alpha")),
    _value1(getParam<Real>("value1")),
    _value2(getParam<Real>("value2"))
{
}

XFEMTwoSideDirichletBC::~XFEMTwoSideDirichletBC() {}

void
XFEMTwoSideDirichletBC::reinitConstraintQuadrature(const ElementPairInfo & element_pair_info)
{
  _interface_normal = element_pair_info._elem1_normal;
  ElemElemConstraint::reinitConstraintQuadrature(element_pair_info);
}

Real
XFEMTwoSideDirichletBC::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0;

  switch (type)
  {
    case Moose::Element:
      r += _alpha * (_u[_qp] - _value1) * _test[_i][_qp];
      break;

    case Moose::Neighbor:
      r += _alpha * (_u_neighbor[_qp] - _value2) * _test_neighbor[_i][_qp];
      break;
  }
  return r;
}

Real
XFEMTwoSideDirichletBC::computeQpJacobian(Moose::DGJacobianType type)
{
  Real r = 0;

  switch (type)
  {
    case Moose::ElementElement:
      r += _alpha * _phi[_j][_qp] * _test[_i][_qp];
      break;

    case Moose::NeighborNeighbor:
      r += _alpha * _phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp];
      break;
  }

  return r;
}
