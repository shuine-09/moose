//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneralPostprocessor.h"

class AlphaThicknessZr;

template <>
InputParameters validParams<AlphaThicknessZr>();

/**
 * Calculates alpha layer thickness given by the C4 model from the oxide/alpha and alpha/beta interfaces positions
 */
class AlphaThicknessZr : public GeneralPostprocessor
{
public:
  static InputParameters validParams();

  AlphaThicknessZr(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual Real getValue() override;

protected:
  /// The alpha layer thickness [um]
  Real _d_alpha;

  /// The oxide/alpha interface position
  const PostprocessorValue & _x_ox_a;

  /// The alpha/beta interface position
  const PostprocessorValue & _x_a_b;
};
