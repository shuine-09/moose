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

#include "DGCHZeroFlux.h"
#include "Function.h"
#include "MooseVariable.h"
#include "libmesh/numeric_vector.h"

#include <cmath>

registerMooseObject("MooseApp", DGCHZeroFlux);

InputParameters
DGCHZeroFlux::validParams()
{
  InputParameters params = IntegratedBC::validParams();
  params.addParam<Real>("alpha",100, "The statblization parameter in the Boundary");
  return params;
}

DGCHZeroFlux::DGCHZeroFlux(const InputParameters & parameters) :
    IntegratedBC(parameters),
    _second_phi(secondPhiFace()),
    _second_test(secondTestFace()),
    _second_u(second()),
    _alpha(getParam<Real>("alpha"))
{
}

Real
DGCHZeroFlux::computeQpResidual()
{
  const unsigned int elem_b_order = static_cast<unsigned int> (_var.order());
  const double h_elem = _current_elem->volume()/_current_side_elem->volume() * 1./std::pow(elem_b_order, 2.);

  Real r = 0;
  r -= (_second_u[_qp].tr() * _grad_test[_i][_qp] * _normals[_qp]);
  r -= (_grad_u[_qp] * _normals[_qp] * _second_test[_i][_qp].tr());
  r += _alpha/h_elem * (_grad_u[_qp] * _normals[_qp] * _grad_test[_i][_qp] * _normals[_qp]);

  return r;
}

Real
DGCHZeroFlux::computeQpJacobian()
{
  const unsigned int elem_b_order = static_cast<unsigned int> (_var.order());
  const double h_elem = _current_elem->volume()/_current_side_elem->volume() * 1./std::pow(elem_b_order, 2.);

  Real r = 0;
  r -= (_second_phi[_j][_qp].tr() * _grad_test[_i][_qp] * _normals[_qp]);
  r -= (_grad_phi[_j][_qp] * _normals[_qp] * _second_test[_i][_qp].tr());
  r += _alpha/h_elem * _grad_phi[_j][_qp] * _normals[_qp] *  _grad_test[_i][_qp] * _normals[_qp];

  return r;
}
