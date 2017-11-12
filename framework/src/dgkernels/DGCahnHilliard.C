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

#include "DGCahnHilliard.h"
#include "NeighborMooseVariableInterface.h"
#include "MooseVariable.h"
#include <cmath>

template <>
InputParameters
validParams<DGCahnHilliard>()
{
  InputParameters params = validParams<DGKernel>();
  params.addRequiredParam<MaterialPropertyName>("kappa_name", "The kappa used with the kernel");
  params.addParam<Real>("eta", 100, "The stablization parameter in DG kernel");
  return params;
}

DGCahnHilliard::DGCahnHilliard(const InputParameters & parameters)
  : DGKernel(parameters),
    _kappa(getMaterialProperty<Real>("kappa_name")),
    _eta(getParam<Real>("eta")),
    _second_phi(secondPhiFace()),
    _second_test(secondTestFace()),
    _second_u(second()),
    _second_phi_neighbor(neighborSecondPhi()),
    _second_test_neighbor(neighborSecondTest()),
    _second_u_neighbor(neighborSecond())
{
}

Real
DGCahnHilliard::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0;

  const unsigned int elem_b_order = static_cast<unsigned int>(_var.order());
  const double h_elem =
      _current_elem->volume() / _current_side_elem->volume() * 1. / std::pow(elem_b_order, 2.);

  switch (type)
  {
    case Moose::Element:
      r += -_kappa[_qp] * 0.5 * (_second_u[_qp].tr() + _second_u_neighbor[_qp].tr()) *
           (_grad_test[_i][_qp] * _normals[_qp]);
      r += -_kappa[_qp] * (_grad_u[_qp] * _normals[_qp] - _grad_u_neighbor[_qp] * _normals[_qp]) *
           0.5 * (_second_test[_i][_qp].tr());
      r += _eta / h_elem * (_grad_u[_qp] * _normals[_qp] - _grad_u_neighbor[_qp] * _normals[_qp]) *
           (_grad_test[_i][_qp] * _normals[_qp]);
      break;

    case Moose::Neighbor:
      r += -_kappa[_qp] * 0.5 * (_second_u[_qp].tr() + _second_u_neighbor[_qp].tr()) *
           (-_grad_test_neighbor[_i][_qp] * _normals[_qp]);
      r += -_kappa[_qp] * (_grad_u[_qp] * _normals[_qp] - _grad_u_neighbor[_qp] * _normals[_qp]) *
           0.5 * (_second_test_neighbor[_i][_qp].tr());
      r += _eta / h_elem * (_grad_u[_qp] * _normals[_qp] - _grad_u_neighbor[_qp] * _normals[_qp]) *
           (-_grad_test_neighbor[_i][_qp] * _normals[_qp]);
      break;
  }

  return r;
}

Real
DGCahnHilliard::computeQpJacobian(Moose::DGJacobianType type)
{
  Real r = 0;

  const unsigned int elem_b_order = static_cast<unsigned int>(_var.order());
  const double h_elem =
      _current_elem->volume() / _current_side_elem->volume() * 1. / std::pow(elem_b_order, 2.);

  switch (type)
  {

    case Moose::ElementElement:
      r += -_kappa[_qp] * 0.5 * (_second_phi[_j][_qp].tr() * _grad_test[_i][_qp] * _normals[_qp]);
      r += -_kappa[_qp] * _grad_phi[_j][_qp] * _normals[_qp] * 0.5 * _second_test[_i][_qp].tr();
      r += _eta / h_elem * (_grad_phi[_j][_qp] * _normals[_qp]) *
           (_grad_test[_i][_qp] * _normals[_qp]);
      break;

    case Moose::ElementNeighbor:
      r += -_kappa[_qp] * 0.5 *
           (_second_phi_neighbor[_j][_qp].tr() * _grad_test[_i][_qp] * _normals[_qp]);
      r += -_kappa[_qp] * (-_grad_phi_neighbor[_j][_qp]) * _normals[_qp] * 0.5 *
           _second_test[_i][_qp].tr();
      r += _eta / h_elem * (-_grad_phi_neighbor[_j][_qp] * _normals[_qp]) *
           (_grad_test[_i][_qp] * _normals[_qp]);
      break;

    case Moose::NeighborElement:
      r += -_kappa[_qp] * 0.5 *
           (_second_phi[_j][_qp].tr() * (-_grad_test_neighbor[_i][_qp]) * _normals[_qp]);
      r += -_kappa[_qp] * _grad_phi[_j][_qp] * _normals[_qp] * 0.5 *
           _second_test_neighbor[_i][_qp].tr();
      r += _eta / h_elem * (_grad_phi[_j][_qp] * _normals[_qp]) *
           (-_grad_test_neighbor[_i][_qp] * _normals[_qp]);
      break;

    case Moose::NeighborNeighbor:
      r += -_kappa[_qp] * 0.5 *
           (_second_phi_neighbor[_j][_qp].tr() * (-_grad_test_neighbor[_i][_qp]) * _normals[_qp]);
      r += -_kappa[_qp] * (-_grad_phi_neighbor[_j][_qp]) * _normals[_qp] * 0.5 *
           _second_test_neighbor[_i][_qp].tr();
      r += _eta / h_elem * (-_grad_phi_neighbor[_j][_qp] * _normals[_qp]) *
           (-_grad_test_neighbor[_i][_qp] * _normals[_qp]);
      break;
  }

  return r;
}
