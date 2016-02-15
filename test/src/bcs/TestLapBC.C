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

#include <cmath>
#include "Assembly.h"

template<>
InputParameters validParams<TestLapBC>()
{
  InputParameters params = validParams<IntegratedBC>();
  return params;
}

TestLapBC::TestLapBC(const InputParameters & parameters) :
    IntegratedBC(parameters),
    _second_phi(_assembly.secondPhiFace()),
    _second_test(_var.secondPhiFace()),
    _second_u(_is_implicit ?_var.secondSln() : _var.secondSlnOld())
{
}

Real
TestLapBC::computeQpResidual()
{
  Real r = 0;
  r += _second_u[_qp].tr() * _test[_i][_qp];
  return r;
}

Real
TestLapBC::computeQpJacobian()
{
  Real r = 0;
  r += _second_phi[_j][_qp].tr() * _test[_i][_qp];
  return r;
}

