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
class INSTwoPhaseSurfaceTension;

declareADValidParams(INSTwoPhaseSurfaceTension);

/**
 *
 */
template <ComputeStage compute_stage>
class INSTwoPhaseSurfaceTension : public ADVectorKernelValue<compute_stage>
{
public:
  INSTwoPhaseSurfaceTension(const InputParameters & parameters);

protected:
  virtual ADRealVectorValue precomputeQpResidual() override;

  const ADVectorVariableValue & _grad_c;
  const ADVariableValue & _curvature;
  const Real _sigma;

  usingVectorKernelValueMembers;
};
