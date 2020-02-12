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

// Forward Declarations
template <ComputeStage>
class DEDMomentSource;

declareADValidParams(DEDMomentSource);

/**
 *
 */
template <ComputeStage compute_stage>
class DEDMomentSource : public ADVectorKernelValue<compute_stage>
{
public:
  DEDMomentSource(const InputParameters & parameters);

protected:
  virtual ADRealVectorValue precomputeQpResidual() override;

  const ADMaterialProperty(RealVectorValue) & _ded_momentum;

  usingVectorKernelValueMembers;
};
