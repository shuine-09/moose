//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PFFractureLineIC.h"
#include "MooseMesh.h"

registerMooseObject("PhaseFieldApp", PFFractureLineIC);

template <>
InputParameters
validParams<PFFractureLineIC>()
{
  InputParameters params = validParams<InitialCondition>();
  params.addRequiredParam<Real>("l",
                                "The length scale parameter used in phase field fracture model.");
  params.addRequiredParam<Point>("point_1", "The coordinate of the first point.");
  params.addRequiredParam<Point>("point_2", "The coordinate of the second point.");
  return params;
}

PFFractureLineIC::PFFractureLineIC(const InputParameters & parameters)
  : InitialCondition(parameters),
    _mesh(_fe_problem.mesh()),
    _l(getParam<Real>("l")),
    _point_1(getParam<Point>("point_1")),
    _point_2(getParam<Point>("point_2"))
{
  if (_mesh.dimension() != 2)
    mooseError("PFFractureLineIC requires a two-dimensional mesh.");
}

Real
PFFractureLineIC::value(const Point & p)
{
  Real min_dist = std::numeric_limits<Real>::max();

  Point c = p - _point_1;
  Point v = (_point_2 - _point_1) / (_point_2 - _point_1).norm();
  Real d = (_point_2 - _point_1).norm();
  Real t = v * c;

  Real dist;
  Point nearest_point;

  if (t < 0)
  {
    dist = (p - _point_1).norm();
    nearest_point = _point_1;
  }
  else if (t > d)
  {
    dist = (p - _point_2).norm();
    nearest_point = _point_2;
  }
  else
  {
    v *= t;
    dist = (p - _point_1 - v).norm();
    nearest_point = (_point_1 + v);
  }

  Point p_nearest_point = nearest_point - p;

  Point normal_ab = Point(-(_point_2 - _point_1)(1), (_point_2 - _point_1)(0), 0);

  if (normal_ab * p_nearest_point < 0)
    dist = -dist;

  if (std::abs(dist) < std::abs(min_dist))
    min_dist = dist;

  return std::exp(-std::abs(min_dist) / _l);
}

RealGradient
PFFractureLineIC::gradient(const Point & /*p*/)
{
  RealGradient gradient = 0.0;
  return gradient;
}
