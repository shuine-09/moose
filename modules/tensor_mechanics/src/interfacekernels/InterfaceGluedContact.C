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

#include "InterfaceGluedContact.h"

#include <cmath>

registerMooseObject("TensorMechanicsApp", InterfaceGluedContact);

template <>
InputParameters
validParams<InterfaceGluedContact>()
{
  InputParameters params = validParams<InterfaceKernel>();
  params.addParam<Real>("alpha", 1e4, "penalty.");
  params.addRequiredParam<unsigned int>("component",
                                        "An integer corresponding to the direction "
                                        "the variable this kernel acts in. (0 for x, "
                                        "1 for y, 2 for z)");
  params.addParam<bool>("tangential", false, "Enforce tangential only.");
  return params;
}

InterfaceGluedContact::InterfaceGluedContact(const InputParameters & parameters)
  : InterfaceKernel(parameters),
    _alpha(getParam<Real>("alpha")),
    _component(getParam<unsigned int>("component")),
    _tangential(getParam<bool>("tangential"))
{
  if (!parameters.isParamValid("boundary"))
  {
    mooseError(
        "In order to use the InterfaceGluedContact dgkernel, you must specify a boundary where "
        "it will live.");
  }
}

Real
InterfaceGluedContact::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0;

  switch (type)
  {
    case Moose::Element:
      r += _alpha * (_u[_qp] - _neighbor_value[_qp]) * _test[_i][_qp];
      break;

    case Moose::Neighbor:
      r -= _alpha * (_u[_qp] - _neighbor_value[_qp]) * _test_neighbor[_i][_qp];
      break;
  }

  Real traction = r;

  if (_tangential)
  {
    if (_component == 0)
      traction = r - r * _normals[_qp](0);
    else if (_component == 1)
      traction = r - r * _normals[_qp](1);
  }

  return traction;
}

Real
InterfaceGluedContact::computeQpJacobian(Moose::DGJacobianType type)
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

  Real dtraction = r;

  if (_tangential)
  {
    if (_component == 0)
      dtraction = r - r * _normals[_qp](0);
    else if (_component == 1)
      dtraction = r - r * _normals[_qp](1);
  }

  return dtraction;
}
