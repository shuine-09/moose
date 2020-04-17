//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADKernel.h"
#include "MeltPoolLevelSetLocation.h"

template <ComputeStage>
class MeltPoolHeatSource;

declareADValidParams(MeltPoolHeatSource);

template <ComputeStage compute_stage>
class MeltPoolHeatSource : public ADKernel<compute_stage>
{
public:
  static InputParameters validParams();

  MeltPoolHeatSource(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual();

  /// Gradient o level set variable
  const ADVectorVariableValue & _grad_ls;

  /// Laser power
  const Real & _power;

  /// Absorption coefficient
  const Real & _alpha;

  /// Laser beam radius
  const Real & _Rb;

  /// Heat transfer coefficient
  const Real & _Ah;

  /// Stefan Boltzmann constant
  const Real & _stefan_boltzmann;

  /// Material emissivity
  const Real & _varepsilon;

  /// Ambient temperature
  const Real & _T0;

  /// Laser position
  const MeltPoolLevelSetLocation & _location;

  usingKernelMembers;
};
