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

#include "DGVesicleShapeDeformation.h"
#include "NeighborMooseVariableInterface.h"
#include "MooseVariable.h"
#include <cmath>

registerMooseObject("MooseApp", DGVesicleShapeDeformation);

template <>
InputParameters
validParams<DGVesicleShapeDeformation>()
{
  InputParameters params = validParams<DGKernel>();
  params.addParam<Real>("eta", 100, "The stablization parameter in DG kernel");
  params.addParam<Real>("epsilon", 0.01, "The interfacial penalty parameter.");
  params.addParam<bool>("rz", false, "The RZ coordindate system.");
  params.addParam<Real>("h_elem", -1.0, "The element size h.");
  return params;
}

DGVesicleShapeDeformation::DGVesicleShapeDeformation(const InputParameters & parameters)
  : DGKernel(parameters),
    _eta(getParam<Real>("eta")),
    _epsilon(getParam<Real>("epsilon")),
    _h_elem(getParam<Real>("h_elem")),
    _rz(getParam<bool>("rz")),
    _second_phi(secondPhiFace()),
    _second_test(secondTestFace()),
    _second_u(second()),
    _second_phi_neighbor(neighborSecondPhi()),
    _second_test_neighbor(neighborSecondTest()),
    _second_u_neighbor(neighborSecond())
{
}

Real
DGVesicleShapeDeformation::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0;

  double h_elem = 0.0;
  if (_h_elem < 0.0)
  {
    const unsigned int elem_b_order = static_cast<unsigned int>(_var.order());
    h_elem =
        _current_elem->volume() / _current_side_elem->volume() * 1. / std::pow(elem_b_order, 2.);
  }
  else
    h_elem = _h_elem;

  Real rz_coord = _q_point[_qp](0);

  Real lap_u = _second_u[_qp].tr();
  Real lap_u_neighbor = _second_u_neighbor[_qp].tr();
  Real lap_test_i = _second_test[_i][_qp].tr();
  Real lap_test_i_neighbor = _second_test_neighbor[_i][_qp].tr();

  if (_rz)
  {
    lap_u += _grad_u[_qp](0) / rz_coord;
    lap_u_neighbor += _grad_u_neighbor[_qp](0) / rz_coord;
    lap_test_i += _grad_test[_i][_qp](0) / rz_coord;
    lap_test_i_neighbor += _grad_test_neighbor[_i][_qp](0) / rz_coord;
  }

  switch (type)
  {
    case Moose::Element:
      r += -0.5 * (lap_u + lap_u_neighbor) * (_grad_test[_i][_qp] * _normals[_qp]);
      r += -(_grad_u[_qp] * _normals[_qp] - _grad_u_neighbor[_qp] * _normals[_qp]) * 0.5 *
           lap_test_i;
      r += _eta / h_elem / _epsilon *
           (_grad_u[_qp] * _normals[_qp] - _grad_u_neighbor[_qp] * _normals[_qp]) *
           (_grad_test[_i][_qp] * _normals[_qp]);
      r *= _epsilon;
      break;

    case Moose::Neighbor:
      r += -0.5 * (lap_u + lap_u_neighbor) * (-_grad_test_neighbor[_i][_qp] * _normals[_qp]);
      r += -(_grad_u[_qp] * _normals[_qp] - _grad_u_neighbor[_qp] * _normals[_qp]) * 0.5 *
           lap_test_i_neighbor;
      r += _eta / h_elem / _epsilon *
           (_grad_u[_qp] * _normals[_qp] - _grad_u_neighbor[_qp] * _normals[_qp]) *
           (-_grad_test_neighbor[_i][_qp] * _normals[_qp]);
      r *= _epsilon;
      break;
  }

  return r;
}

Real
DGVesicleShapeDeformation::computeQpJacobian(Moose::DGJacobianType type)
{
  Real r = 0;

  double h_elem = 0.0;
  if (_h_elem < 0.0)
  {
    const unsigned int elem_b_order = static_cast<unsigned int>(_var.order());
    h_elem =
        _current_elem->volume() / _current_side_elem->volume() * 1. / std::pow(elem_b_order, 2.);
  }
  else
    h_elem = _h_elem;

  Real rz_coord = _q_point[_qp](0);

  Real lap_test_i = _second_test[_i][_qp].tr();
  Real lap_test_i_neighbor = _second_test_neighbor[_i][_qp].tr();
  Real lap_phi_j = _second_phi[_j][_qp].tr();
  Real lap_phi_j_neighbor = _second_phi_neighbor[_j][_qp].tr();

  if (_rz)
  {
    lap_test_i += _grad_test[_i][_qp](0) / rz_coord;
    lap_test_i_neighbor += _grad_test_neighbor[_i][_qp](0) / rz_coord;
    lap_phi_j += _grad_phi[_j][_qp](0) / rz_coord;
    lap_phi_j_neighbor += _grad_phi_neighbor[_j][_qp](0) / rz_coord;
  }

  switch (type)
  {
    case Moose::ElementElement:
      r += -0.5 * (lap_phi_j * _grad_test[_i][_qp] * _normals[_qp]);
      r += -_grad_phi[_j][_qp] * _normals[_qp] * 0.5 * lap_test_i;
      r += _eta / h_elem / _epsilon * (_grad_phi[_j][_qp] * _normals[_qp]) *
           (_grad_test[_i][_qp] * _normals[_qp]);
      r *= _epsilon;
      break;

    case Moose::ElementNeighbor:
      r += -0.5 * (lap_phi_j_neighbor * _grad_test[_i][_qp] * _normals[_qp]);
      r += -(-_grad_phi_neighbor[_j][_qp]) * _normals[_qp] * 0.5 * lap_test_i;
      r += _eta / h_elem / _epsilon * (-_grad_phi_neighbor[_j][_qp] * _normals[_qp]) *
           (_grad_test[_i][_qp] * _normals[_qp]);
      r *= _epsilon;
      break;

    case Moose::NeighborElement:
      r += -0.5 * (lap_phi_j * (-_grad_test_neighbor[_i][_qp]) * _normals[_qp]);
      r += -_grad_phi[_j][_qp] * _normals[_qp] * 0.5 * lap_test_i_neighbor;
      r += _eta / h_elem / _epsilon * (_grad_phi[_j][_qp] * _normals[_qp]) *
           (-_grad_test_neighbor[_i][_qp] * _normals[_qp]);
      r *= _epsilon;
      break;

    case Moose::NeighborNeighbor:
      r += -0.5 * (lap_phi_j_neighbor * (-_grad_test_neighbor[_i][_qp]) * _normals[_qp]);
      r += -(-_grad_phi_neighbor[_j][_qp]) * _normals[_qp] * 0.5 * lap_test_i_neighbor;
      r += _eta / h_elem / _epsilon * (-_grad_phi_neighbor[_j][_qp] * _normals[_qp]) *
           (-_grad_test_neighbor[_i][_qp] * _normals[_qp]);
      r *= _epsilon;
      break;
  }

  return r;
}
