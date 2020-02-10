//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADVectorKernel.h"

template <ComputeStage>
class LevelSetGradientRegulization;

declareADValidParams(LevelSetGradientRegulization);

template <ComputeStage compute_stage>
class LevelSetGradientRegulization : public ADVectorKernel<compute_stage>
{
public:
  static InputParameters validParams();

  LevelSetGradientRegulization(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const ADVariableGradient & _grad_c;

  const bool _normal_regulization;

  usingVectorKernelMembers;
};
