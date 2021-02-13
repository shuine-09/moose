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

class WeightGainSpaceIntegralZr;

template <>
InputParameters validParams<WeightGainSpaceIntegralZr>();

/**
 * Calculates weight gain given by the C4 model from the space integral of the oxygen concentration
 */
class WeightGainSpaceIntegralZr : public GeneralPostprocessor
{
public:
  static InputParameters validParams();

  WeightGainSpaceIntegralZr(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual Real getValue() override;

protected:
  /// Temperature [K] to get the interfaces concentrations
  Real _temperature;
  /// The weight gain [mg/cm^2]
  Real _wg;

  /// The oxygen concentration (weak and reduced) integral [um² because integrated over x AND y]
  const PostprocessorValue & _C_integral;

  /// The y dimension of the mesh  [um]
  Real _ymax;

  /// The oxide thickness [um]
  const PostprocessorValue & _delta;

  /// The alpha layer thickness [um]
  const PostprocessorValue & _d_alpha;
};
