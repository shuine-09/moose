//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// MOOSE includes
#include "LevelSetOlssonBubble.h"

registerMooseObject("LevelSetApp", LevelSetOlssonBubble);

template <>
InputParameters
validParams<LevelSetOlssonBubble>()
{
  InputParameters params = validParams<Function>();
  params.addClassDescription("Implementation of 'bubble' ranging from 0 to 1.");
  params.addParam<RealVectorValue>(
      "center", RealVectorValue(0.5, 0.5, 0), "The center of the bubble.");
  params.addParam<Real>("radius", 0.15, "The radius of the bubble.");
  params.addParam<Real>("epsilon", 0.01, "The interface thickness.");
  return params;
}

LevelSetOlssonBubble::LevelSetOlssonBubble(const InputParameters & parameters)
  : Function(parameters),
    _center(getParam<RealVectorValue>("center")),
    _radius(getParam<Real>("radius")),
    _epsilon(getParam<Real>("epsilon"))
{
}

Real
LevelSetOlssonBubble::value(Real /*t*/, const Point & z) const
{
  // Real _a = 0.3125;
  // Real _b = 0.2;
  //
  // Point p(std::abs(z(0)), std::abs(z(1)), std::abs(z(2)));
  // Real a = _a;
  // Real b = _b;
  // if (p(0) > p(1))
  // {
  //   Real temp = p(1);
  //   p(1) = p(0);
  //   p(0) = temp;
  //   a = _b;
  //   b = _a;
  // }
  //
  // Real l = b * b - a * a;
  // Real m = a * p(0) / l;
  // Real m2 = m * m;
  // Real n = b * p(1) / l;
  // Real n2 = n * n;
  // Real c = (m2 + n2 - 1.0) / 3.0;
  // Real c3 = c * c * c;
  // Real q = c3 + m2 * n2 * 2.0;
  // Real d = c3 + m2 * n2;
  // Real g = m + m * n2;
  //
  // Real co;
  //
  // if (d < 0.0)
  // {
  //   Real p = std::acos(q / c3) / 3.0;
  //   Real s = std::cos(p);
  //   Real t = std::sin(p) * sqrt(3.0);
  //   Real rx = std::sqrt(-c * (s + t + 2.0) + m2);
  //   Real ry = std::sqrt(-c * (s - t + 2.0) + m2);
  //   co = (ry + std::copysignf(1.0, l) * rx + std::abs(g) / (rx * ry) - m) / 2.0;
  // }
  // else
  // {
  //   Real h = 2.0 * m * n * std::sqrt(d);
  //   Real s = std::copysignf(1.0, q + h) * std::pow(std::abs(q + h), 1.0 / 3.0);
  //   Real u = std::copysignf(1.0, q - h) * std::pow(std::abs(q - h), 1.0 / 3.0);
  //   Real rx = -s - u - c * 4.0 + 2.0 * m2;
  //   Real ry = (s - u) * std::sqrt(3.0);
  //   Real rm = std::sqrt(rx * rx + ry * ry);
  //   Real p = ry / std::sqrt(rm - rx);
  //   co = (p + 2.0 * g / rm - m) / 2.0;
  // }
  //
  // Real si = std::sqrt(1.0 - co * co);
  //
  // Point closestPoint(a * co, b * si, 0);
  //
  // Real dist = (closestPoint - p).size() * std::copysignf(1.0, p(1) - closestPoint(1));
  //
  // const Real x = dist / _epsilon;
  // return 1.0 / (1 + std::exp(x));

  // const Real x = ((z - _center).norm() - _radius) / _epsilon;
  // return 1.0 / (1 + std::exp(x));

  const Real dist1 = (z - _center).norm() - _radius;
  const Real dist2 = z(1) - 2.0;

  const Real dist = std::min(dist1, dist2);
  const Real x = dist / _epsilon;
  return 1.0 / (1 + std::exp(x));
}

RealGradient
LevelSetOlssonBubble::gradient(Real /*t*/, const Point & p) const
{
  Real norm = (p - _center).norm();
  Real g = (norm - _radius) / _epsilon;
  RealGradient output;

  Real g_prime;
  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  {
    g_prime = (p(i) - _center(i)) / (_epsilon * norm);
    output(i) = (g_prime * std::exp(g)) / ((std::exp(g) + 1) * (std::exp(g) + 1));
  }
  return output;
}
