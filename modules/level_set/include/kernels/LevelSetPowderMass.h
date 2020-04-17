//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADKernelValue.h"
#include "MeltPoolLevelSetLocation.h"

// Forward declarations
template <ComputeStage>
class LevelSetPowderMass;

declareADValidParams(LevelSetPowderMass);

/**
 */

template <ComputeStage compute_stage>
class LevelSetPowderMass : public ADKernelValue<compute_stage>
{
public:
  static InputParameters validParams();

  LevelSetPowderMass(const InputParameters & parameters);

protected:
  virtual ADReal precomputeQpResidual() override;

  /// Density
  const ADMaterialProperty(Real) & _rho;

  const ADVectorVariableValue & _grad_ls;

  /// Mass rate
  const Real & _mass_rate;

  /// Mass radius
  const Real & _mass_radius;

  /// Laser position
  const MeltPoolLevelSetLocation & _location;

  usingKernelValueMembers;
};
