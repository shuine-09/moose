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

registerMooseObject("MooseApp", LinearStommelMunk);

InputParameters
LinearStommelMunk::validParams()
{
  InputParameters p = Kernel::validParams();
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

Real
LinearStommelMunk::computeQpResidual()
{
  return _eps_s * _grad_u[_qp] * _grad_test[_i][_qp] + _eps_m * _second_u[_qp].tr() * _second_test[_i][_qp].tr() - _grad_u[_qp](0) * _test[_i][_qp];
  //std::cout << "qp = " << _qp << ", point = " << _q_point[_qp] << ", i = " << _i << std::endl;
  //std::cout << "i = " << _i << ", test_function = " << _test[_i][_qp] << std::endl;
  //std::cout << "gradient = " << _grad_test[_i][_qp] << std::endl;
  //std::cout << _second_test[_i][_qp] << std::endl;

  //return _second_u[_qp].tr() * _second_test[_i][_qp].tr();

  //return 0.0;
}

Real
LinearStommelMunk::computeQpJacobian()
{
  return _eps_s * _grad_phi[_j][_qp] * _grad_test[_i][_qp] + _eps_m * _second_phi[_j][_qp].tr() * _second_test[_i][_qp].tr() - _grad_phi[_j][_qp](0) * _test[_i][_qp];
  //return  _second_phi[_j][_qp].tr() * _second_test[_i][_qp].tr();
  //return 0.0;
}
