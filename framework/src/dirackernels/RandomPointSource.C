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

#include "RandomPointSource.h"
#include "MooseRandom.h"
#include "MooseMesh.h"

template<>
InputParameters validParams<RandomPointSource>()
{
  InputParameters params = validParams<DiracKernel>();
  params.addRequiredParam<Real>("value", "The value of the point source");
  params.addParam<unsigned int>("num_srcs", 1, "the number of point sources each time step");
  params.addParam<Real>("dt_one", 1, "Time step size for one hit");

  params.addRequiredParam<Real>("x1", "The x coordinate of the lower left-hand corner of the box");
  params.addRequiredParam<Real>("y1", "The y coordinate of the lower left-hand corner of the box");
  params.addParam<Real>("z1", 0.0, "The z coordinate of the lower left-hand corner of the box");

  params.addRequiredParam<Real>("x2", "The x coordinate of the upper right-hand corner of the box");
  params.addRequiredParam<Real>("y2", "The y coordinate of the upper right-hand corner of the box");
  params.addParam<Real>("z2", 0.0, "The z coordinate of the upper right-hand corner of the box");

  return params;
}

RandomPointSource::RandomPointSource(const InputParameters & parameters) :
    DiracKernel(parameters),
    _mesh_dimension(_mesh.dimension()),
    _value(getParam<Real>("value")),
    _num_srcs(getParam<unsigned int>("num_srcs")),
    _dt_one(getParam<Real>("dt_one")),
    _x1(getParam<Real>("x1")),
    _y1(getParam<Real>("y1")),
    _z1(getParam<Real>("z1")),
    _x2(getParam<Real>("x2")),
    _y2(getParam<Real>("y2")),
    _z2(getParam<Real>("z2")),
    _bottom_left(_x1,_y1,_z1),
    _top_right(_x2,_y2,_z2),
    _range(_top_right - _bottom_left)
{
}

void
RandomPointSource::addPoints()
{
  MooseRandom::seed(_t_step + 1111);

  unsigned int num_dirac_srcs = _num_srcs;

  // Change with time step size
  if (_num_srcs == 0)
  {
    num_dirac_srcs = _dt / _dt_one + 0.5;
    if (num_dirac_srcs < 1)
      num_dirac_srcs = 1;
  }

  for (unsigned int i = 0; i < num_dirac_srcs; ++i)
  {
    _p = 0.0;

    for (unsigned int j = 0; j < 2; ++j)
    {
      if (j<_mesh_dimension)
        _p(j) = _bottom_left(j) + _range(j) * MooseRandom::rand();
    }

    addPoint(_p);
  }
}

Real
RandomPointSource::computeQpResidual()
{
  // This is negative because it's a forcing function that has been brought over to the left side
  return -_test[_i][_qp] * _value;
}
