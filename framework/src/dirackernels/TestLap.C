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

#include "TestLap.h"

registerMooseObject("MooseApp", TestLap);

InputParameters
TestLap::validParams()
{
  InputParameters params = DiracKernel::validParams();
  return params;
}

TestLap::TestLap(const InputParameters & parameters) :
    DiracKernel(parameters),
    _second_phi(secondPhi()),
    _second_test(secondTest()),
    _second_u(second())
{
}

void
TestLap::addPoints()
{
  Point p(1.0, 0.5, 0);
  addPoint(_current_elem,p);
}

Real
TestLap::computeQpResidual()
{
//  This is negative because it's a forcing function that has been brought over to the left side
  std::cout << "TESTLAP::point = " << _current_point << ", i = " << _i << std::endl;
  std::cout << "grad = " << _grad_test[_i][_qp] << std::endl;
  std::cout << _second_test[_i][_qp] << std::endl;
  return 0.0;
  //return -_test[_i][_qp]*_value;
}
