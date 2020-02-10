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
#include "ADKernel.h"

// Forward declarations
template <ComputeStage>
class LevelSetCurvatureRegulization;

declareADValidParams(LevelSetCurvatureRegulization);

/**
 *
 */
template <ComputeStage compute_stage>
class LevelSetCurvatureRegulization : public ADKernel<compute_stage>
{
public:
  LevelSetCurvatureRegulization(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual();

  const ADVectorVariableValue & _grad_c;

  const Real _epsilon;

  usingKernelMembers;
};
