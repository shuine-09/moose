//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADKernelGrad.h"

template <ComputeStage>
class HeatConvection;

declareADValidParams(HeatConvection);

template <ComputeStage compute_stage>
class HeatConvection : public ADKernelGrad<compute_stage>
{
public:
  static InputParameters validParams();

  HeatConvection(const InputParameters & parameters);

protected:
  virtual ADRealVectorValue precomputeQpResidual() override;

  const ADVectorVariableValue & _velocity;
  const ADMaterialProperty(Real) & _rho;
  const ADMaterialProperty(Real) & _h;
  const ADVariableValue & _ls;

  usingKernelGradMembers;
};
