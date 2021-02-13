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

class WeightGainZr;

template <>
InputParameters validParams<WeightGainZr>();

/**
 * Calculates weight gain given by the C4 model from the time integral of the vacancy flux
 */
class WeightGainZr : public GeneralPostprocessor
{
public:
  static InputParameters validParams();

  WeightGainZr(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual Real getValue() override;

protected:
  /// The weight gain [mg/cm^2]
  Real _wg;

  /// Temperature [K] to get the "initial" weight gain at 20s
  Real _temperature;

  /// The vacancy flux integral [/m^2]
  const PostprocessorValue & _flux_integral;
};
