//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE includes
#include "Function.h"

class LevelSetOlssonPlane;

template <>
InputParameters validParams<LevelSetOlssonPlane>();

/**
 * Implements the "bubble" function from Olsson and Kreiss (2005).
 */
class LevelSetOlssonPlane : public Function
{
public:
  LevelSetOlssonPlane(const InputParameters & parameters);

  virtual Real value(Real /*t*/, const Point & p) const override;

  virtual RealGradient gradient(Real /*t*/, const Point & p) const override;

protected:
  /// The interface thickness
  const Real & _epsilon;
  const bool & _z;
};
