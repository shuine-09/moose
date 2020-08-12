//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "InitialCondition.h"

class PhaseFieldDamageIC;

class PhaseFieldDamageIC : public InitialCondition
{
public:
  static InputParameters validParams();
  PhaseFieldDamageIC(const InputParameters & parameters);

protected:
  virtual Real dist(const Point & p);
  virtual Real value(const Point & p);

  std::vector<Real> _x1;
  std::vector<Real> _y1;
  std::vector<Real> _z1;
  std::vector<Real> _x2;
  std::vector<Real> _y2;
  std::vector<Real> _z2;
  const Real _d0;
  const Real _l;
};
