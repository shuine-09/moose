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

#include "DGLinearStommelMunk.h"
#include "NeighborMooseVariableInterface.h"
#include "MooseVariable.h"
#include <cmath>

registerMooseObject("MooseApp", DGLinearStommelMunk);

InputParameters
DGLinearStommelMunk::validParams()
{
  InputParameters params = DGKernel::validParams();
  params.addParam<Real>("eta", 100, "The stablization parameter in DG kernel");
  return params;
}

DGLinearStommelMunk::DGLinearStommelMunk(const InputParameters & parameters) :
    DGKernel(parameters),
    _eta(getParam<Real>("eta")),
    _second_phi(secondPhiFace()),
    _second_test(secondTestFace()),
    _second_u(second()),
    _second_phi_neighbor(neighborSecondPhi()),
    _second_test_neighbor(neighborSecondTest()),
    _second_u_neighbor(neighborSecond())
{
  _eps_s = 0.05;
  _eps_m = 6.0e-5;
}

Real
DGLinearStommelMunk::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0;

  const unsigned int elem_b_order = static_cast<unsigned int> (_var.order());
  const double h_elem = _current_elem->volume()/_current_side_elem->volume() * 1./std::pow(elem_b_order, 2.);

  switch (type)
  {
  case Moose::Element:
    r += -_eps_m * 0.5 * (_second_u[_qp].tr() + _second_u_neighbor[_qp].tr()) * (_grad_test[_i][_qp] * _normals[_qp]);
    r += -_eps_m * (_grad_u[_qp] * _normals[_qp] - _grad_u_neighbor[_qp] * _normals[_qp]) * 0.5 * (_second_test[_i][_qp].tr());
    r += _eps_m * _eta / h_elem * (_grad_u[_qp] * _normals[_qp] - _grad_u_neighbor[_qp] * _normals[_qp]) * (_grad_test[_i][_qp] * _normals[_qp]);
    break;

  case Moose::Neighbor:
    r += -_eps_m * 0.5 * (_second_u[_qp].tr() + _second_u_neighbor[_qp].tr()) * (-_grad_test_neighbor[_i][_qp] * _normals[_qp]);
    r += -_eps_m * (_grad_u[_qp] * _normals[_qp] - _grad_u_neighbor[_qp] * _normals[_qp]) * 0.5 * (_second_test_neighbor[_i][_qp].tr());
    r += _eps_m * _eta / h_elem * (_grad_u[_qp] * _normals[_qp] - _grad_u_neighbor[_qp] * _normals[_qp]) * (-_grad_test_neighbor[_i][_qp] * _normals[_qp]);
    break;
  }

  return r;
}

Real
DGLinearStommelMunk::computeQpJacobian(Moose::DGJacobianType type)
{
  Real r = 0;

  const unsigned int elem_b_order = static_cast<unsigned int> (_var.order());
  const double h_elem = _current_elem->volume()/_current_side_elem->volume() * 1./std::pow(elem_b_order, 2.);

  switch (type)
  {

  case Moose::ElementElement:
    r += -_eps_m * 0.5 * (_second_phi[_j][_qp].tr() * _grad_test[_i][_qp] * _normals[_qp]);
    r += -_eps_m * _grad_phi[_j][_qp] * _normals[_qp] * 0.5 * _second_test[_i][_qp].tr();
    r += _eps_m * _eta / h_elem * (_grad_phi[_j][_qp] * _normals[_qp]) * (_grad_test[_i][_qp] * _normals[_qp]);
    break;

  case Moose::ElementNeighbor:
    r += -_eps_m * 0.5 * (_second_phi_neighbor[_j][_qp].tr() * _grad_test[_i][_qp] * _normals[_qp]);
    r += -_eps_m * (-_grad_phi_neighbor[_j][_qp]) * _normals[_qp] * 0.5 * _second_test[_i][_qp].tr();
    r += _eps_m * _eta / h_elem * (-_grad_phi_neighbor[_j][_qp] * _normals[_qp]) * (_grad_test[_i][_qp] * _normals[_qp]);
    break;

  case Moose::NeighborElement:
    r += -_eps_m * 0.5 * (_second_phi[_j][_qp].tr() * (-_grad_test_neighbor[_i][_qp]) * _normals[_qp]);
    r += -_eps_m * _grad_phi[_j][_qp] * _normals[_qp] * 0.5 * _second_test_neighbor[_i][_qp].tr();
    r += _eps_m * _eta / h_elem * (_grad_phi[_j][_qp] * _normals[_qp]) * (-_grad_test_neighbor[_i][_qp] * _normals[_qp]);
    break;

  case Moose::NeighborNeighbor:
    r += -_eps_m * 0.5 * (_second_phi_neighbor[_j][_qp].tr() * (-_grad_test_neighbor[_i][_qp]) * _normals[_qp]);
    r += -_eps_m * (-_grad_phi_neighbor[_j][_qp]) * _normals[_qp] * 0.5 * _second_test_neighbor[_i][_qp].tr();
    r += _eps_m * _eta / h_elem * (-_grad_phi_neighbor[_j][_qp] * _normals[_qp]) * (-_grad_test_neighbor[_i][_qp] * _normals[_qp]);
    break;
  }

  return r;
}
