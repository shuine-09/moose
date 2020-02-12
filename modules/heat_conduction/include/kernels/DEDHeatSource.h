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

template <ComputeStage>
class DEDHeatSource;
class Function;

declareADValidParams(DEDHeatSource);

template <ComputeStage compute_stage>
class DEDHeatSource : public ADKernel<compute_stage>
{
public:
  static InputParameters validParams();

  DEDHeatSource(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual();

  const ADMaterialProperty(Real) & _heat_source;

  usingKernelMembers;
};
