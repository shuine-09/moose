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

#include "LinearStommelMunk.h"

template<>
InputParameters validParams<LinearStommelMunk>()
{
  InputParameters p = validParams<Kernel>();
  return p;
}

LinearStommelMunk::LinearStommelMunk(const InputParameters & parameters) :
    Kernel(parameters),
    _second_phi(secondPhi()),
    _second_test(secondTest()),
    _second_u(second())
{
  _eps_s = 0.05;
  _eps_m = 6.0e-5;
}

LinearStommelMunk::~LinearStommelMunk()
{
}

Real
LinearStommelMunk::computeQpResidual()
{
  return _eps_s * _grad_u[_qp] * _grad_test[_i][_qp] + _eps_m * _second_u[_qp].tr() * _second_test[_i][_qp].tr() - _grad_u[_qp](0) * _test[_i][_qp];
}

Real
LinearStommelMunk::computeQpJacobian()
{
  return _eps_s * _grad_phi[_j][_qp] * _grad_test[_i][_qp] + _eps_m * _second_phi[_j][_qp].tr() * _second_test[_i][_qp].tr() - _grad_phi[_j][_qp](0) * _test[_i][_qp];
}

