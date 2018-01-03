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

#include "InterfaceCohesiveZone.h"

#include <cmath>

template <>
InputParameters
validParams<InterfaceCohesiveZone>()
{
  InputParameters params = validParams<InterfaceKernel>();
  return params;
}

InterfaceCohesiveZone::InterfaceCohesiveZone(const InputParameters & parameters)
  : InterfaceKernel(parameters)
{
  if (!parameters.isParamValid("boundary"))
  {
    mooseError(
        "In order to use the InterfaceCohesiveZone dgkernel, you must specify a boundary where "
        "it will live.");
  }
  _max_normal_separation = &getMaterialProperty<Real>("max_normal_separation");
  _max_normal_separation_old = &getMaterialPropertyOld<Real>("max_normal_separation");
}

Real
InterfaceCohesiveZone::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0;

  switch (type)
  {
    case Moose::Element:
      r += 1.0e8 * (_u[_qp] - _neighbor_value[_qp]) * _test[_i][_qp];
      break;

    case Moose::Neighbor:
      r -= 1.0e8 * (_u[_qp] - _neighbor_value[_qp]) * _test_neighbor[_i][_qp];
      break;
  }
  return r;
}

Real
InterfaceCohesiveZone::computeQpJacobian(Moose::DGJacobianType type)
{
  Real r = 0;

  switch (type)
  {
    case Moose::ElementElement:
      r += 1.0e8 * _phi[_j][_qp] * _test[_i][_qp];
      break;

    case Moose::ElementNeighbor:
      r -= 1.0e8 * _phi_neighbor[_j][_qp] * _test[_i][_qp];
      break;

    case Moose::NeighborElement:
      r -= 1.0e8 * _phi[_j][_qp] * _test_neighbor[_i][_qp];
      break;

    case Moose::NeighborNeighbor:
      r += 1.0e8 * _phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp];
      break;
  }

  return r;
}
