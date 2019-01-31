//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef ADCOMPUTEISOTROPICLINEARELASTICPFFRACTURESTRESS_H
#define ADCOMPUTEISOTROPICLINEARELASTICPFFRACTURESTRESS_H

#include "ADComputeStressBase.h"
#include "ADMaterial.h"

// Forward Declarations
template <ComputeStage>
class ADComputeIsotropicLinearElasticPFFractureStress;

declareADValidParams(ADComputeIsotropicLinearElasticPFFractureStress);

/**
 * ADComputeIsotropicLinearElasticPFFractureStress computes the stress following linear elasticity
 * theory (small strains)
 */
template <ComputeStage compute_stage>
class ADComputeIsotropicLinearElasticPFFractureStress : public ADComputeStressBase<compute_stage>
{
public:
  ADComputeIsotropicLinearElasticPFFractureStress(const InputParameters & parameters);

protected:
  virtual void computeQpStress();

  /// Coupled order parameter defining the crack
  const ADVariableValue & _c;

  const VariableValue & _c_old;

  /// Small number to avoid non-positive definiteness at or near complete damage
  const Real _kdamage;

  /// Use current value of history variable
  bool _use_current_hist;

  /// Elastic energy and derivatives, declared in this material

  /// History variable that prevents crack healing, declared in this material
  ADMaterialProperty(Real) & _hist;

  /// Old value of history variable
  const MaterialProperty<Real> & _hist_old;

  usingComputeStressBaseMembers;
};

#endif // ADCOMPUTEISOTROPICLINEARELASTICPFFRACTURESTRESS_H
