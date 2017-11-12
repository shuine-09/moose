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

#include "TestLapBC.h"
#include "Function.h"

#include "libmesh/numeric_vector.h"

registerMooseObject("MooseApp", TestLapBC);

#include <cmath>

InputParameters
TestLapBC::validParams()
{
  InputParameters params = IntegratedBC::validParams();

  return params;
}

TestLapBC::TestLapBC(const InputParameters & parameters) :
    IntegratedBC(parameters),
    _second_phi(secondPhi()),
    _second_test(secondTest()),
    _second_u(second())
{
}

Real
TestLapBC::computeQpResidual()
{
  Real r = 0;
  r += _second_u[_qp].tr() * _second_test[_i][_qp].tr();

  //std::cout << "u = " << _u[_qp] << std::endl;
  //std::cout << "grad_u = " << _grad_u[_qp] << std::endl;
  //std::cout << "second_u = " << _second_u[_qp] << std::endl;

  std::cout << "qp = " << _qp << ", point = " << _q_point[_qp] << ", i = " << _i << std::endl;
  //std::cout << "i = " << _i << ", test_function = " << _test[_i][_qp] << std::endl;
  std::cout << "gradient = " << _grad_test[_i][_qp] << std::endl;
  std::cout << _second_test[_i][_qp] << std::endl;

  return r;
}

Real
TestLapBC::computeQpJacobian()
{
  Real r = 0;
  r += _second_phi[_j][_qp].tr() * _second_test[_i][_qp].tr();

  //std::cout << "j = " << _j << ", i = " << _i << std::endl;
  //std::cout << _second_phi[_j][_qp] << std::endl;
  //std::cout << _second_test[_i][_qp] << std::endl;

  //std::cout << "(" << _j << ", " << _i << ") = " << r << std::endl;

  return r;
}
