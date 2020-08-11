//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "DirichletBCBase.h"

class DirichletBCRightC4Zr;

template <>
InputParameters validParams<DirichletBCRightC4Zr>();

/**
 * Boundary condition of a Dirichlet type
 *
 * Sets the value in the node
 */
class DirichletBCRightC4Zr: public DirichletBCBase
{
public:
  static InputParameters validParams();

  DirichletBCRightC4Zr(const InputParameters & parameters);

protected:
  virtual Real computeQpValue() override;

  /// Boolean specifying if there are 2 interfaces (oxide/alpha and alpha/beta) or only one
  bool _two_interfaces;
  /// The value for this BC
  Real _value;

  /// Temperature [K] : fixes the jump values at interface thus the right boundary (weak discontinuity)
  Real _temperature_aox;
  Real _temperature_ab;


};
