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

class OxideThicknessZr;

template <>
InputParameters validParams<OxideThicknessZr>();

/**
 * Calculates oxide thickness given by the C4 model from the oxide/alpha interface position
 */
class OxideThicknessZr : public GeneralPostprocessor
{
public:
  static InputParameters validParams();

  OxideThicknessZr(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual Real getValue() override;

protected:
  /// The oxide thickness [um]
  Real _delta;

  /// The oxide/alpha interface position
  const PostprocessorValue & _x_ox_a;
};
