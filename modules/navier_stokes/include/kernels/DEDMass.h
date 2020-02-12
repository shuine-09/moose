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
class DEDMass;

declareADValidParams(DEDMass);

/**
 * This class computes the mass equation residual and Jacobian
 * contributions (the latter using automatic differentiation) for the incompressible Navier-Stokes
 * equations.
 */
template <ComputeStage compute_stage>
class DEDMass : public ADKernelValue<compute_stage>
{
public:
  DEDMass(const InputParameters & parameters);

protected:
  ADReal precomputeQpResidual() override;

  /// The strong residual of the mass equation, computed using INSADMaterial
  const ADMaterialProperty(Real) & _mass_strong_residual;

  const ADMaterialProperty(Real) & _rho;
  const ADMaterialProperty(RealVectorValue) & _grad_rho;
  const ADMaterialProperty(Real) & _mass_change;
  const ADMaterialProperty(Real) & _powder_feed;
  const ADVectorVariableValue & _grad_ls;
  usingKernelValueMembers;
};
