/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "XFEMEqualValueConstraint.h"
#include "FEProblem.h"
#include "Assembly.h"

// libMesh includes
#include "libmesh/quadrature.h"

template<>
InputParameters validParams<XFEMEqualValueConstraint>()
{
  InputParameters params = validParams<XFEMElementConstraint>();
  return params;
}

XFEMEqualValueConstraint::XFEMEqualValueConstraint(const InputParameters & parameters) :
    XFEMElementConstraint(parameters)
{
}

XFEMEqualValueConstraint::~XFEMEqualValueConstraint()
{
}

Real 
XFEMEqualValueConstraint::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0;
  
  switch (type)
  {
    case Moose::Element:
      r += 100000 * (_u[_qp] - _u_neighbor[_qp]) * _test[_i][_qp];
      break;

    case Moose::Neighbor:
      r -= 100000 * (_u[_qp] - _u_neighbor[_qp]) * _test_neighbor[_i][_qp];
      break;
  }
  return r;
}

Real 
XFEMEqualValueConstraint::computeQpJacobian(Moose::DGJacobianType type)
{
  Real r = 0;

  switch (type)
  {

    case Moose::ElementElement:
     r += 100000 * _phi[_j][_qp] * _test[_i][_qp];
      break;

    case Moose::ElementNeighbor:
     r -= 100000 * _phi_neighbor[_j][_qp] * _test[_i][_qp];
      break;

    case Moose::NeighborElement:
     r -= 100000 * _phi[_j][_qp] * _test_neighbor[_i][_qp];
      break;

    case Moose::NeighborNeighbor:
     r += 100000 * _phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp];
      break;
  }

  return r;
}
