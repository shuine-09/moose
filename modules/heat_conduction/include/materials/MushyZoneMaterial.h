/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#pragma once

#include "ADMaterial.h"

#define usingMushyZoneMateriallMembers usingMaterialMembers;

template <ComputeStage>
class MushyZoneMaterial;

declareADValidParams(MushyZoneMaterial);

/**
 *
 */

template <ComputeStage compute_stage>
class MushyZoneMaterial : public ADMaterial<compute_stage>
{
public:
  MushyZoneMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// Temperature variable
  const ADVariableValue & _temp;

  /// Solidus temperature
  const Real & _solidus_temperature;

  /// Liquidus  temperature
  const Real & _liquidus_temperature;

  /// Liquid mass fraction
  ADMaterialProperty(Real) & _f_l;

  /// Solid mass fraction
  ADMaterialProperty(Real) & _f_s;

  /// Liquid volume fraction
  ADMaterialProperty(Real) & _g_l;

  /// Solid volume fraction
  ADMaterialProperty(Real) & _g_s;

  /// Solid density
  const Real & _rho_s;

  /// Liquid density
  const Real & _rho_l;

  usingMaterialMembers;
};
