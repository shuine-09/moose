//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#ifndef ADPFFRACTURE_H
#define ADPFFRACTURE_H

#include "ADKernel.h"

// Forward Declarations
template <ComputeStage>
class ADPFFracture;

declareADValidParams(ADPFFracture);

/**
 * Phase field based fracture model
 * This kernel computes the residual and jacobian for bulk free energy contribution to c
 * Refer to Formulation: Miehe et. al., Int. J. Num. Methods Engg., 2010, 83. 1273-1311 Equation
 * 63
 */
template <ComputeStage compute_stage>
class ADPFFracture : public ADKernel<compute_stage>
{
public:
  ADPFFracture(const InputParameters & parameters);

protected:
  virtual ADResidual computeQpResidual() override;

  /// Critical energy release rate for fracture
  const MaterialProperty<Real> & _gc_prop;

  /// Characteristic length, controls damage zone thickness
  const MaterialProperty<Real> & _l;

  const ADMaterialProperty(Real) & _hist;

  const MaterialProperty<Real> & _hist_old;

  usingKernelMembers;
  //@}
};
#endif // ADPFFRACTURE_H
