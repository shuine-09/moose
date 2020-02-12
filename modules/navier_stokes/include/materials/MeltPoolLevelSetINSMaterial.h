/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#pragma once

#include "INSADTauMaterial.h"
#include "MeltPoolLevelSetLocation.h"

#define usingMeltPoolLevelSetINSMateriallMembers usingMaterialMembers;

template <ComputeStage>
class MeltPoolLevelSetINSMaterial;

declareADValidParams(MeltPoolLevelSetINSMaterial);

/**
 *
 */
template <ComputeStage compute_stage>
class MeltPoolLevelSetINSMaterial : public INSADTauMaterial<compute_stage>
{
public:
  MeltPoolLevelSetINSMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// Level set variable
  const ADVariableValue & _ls;

  /// Gradient o level set variable
  const ADVariableGradient & _grad_ls;

  /// Temperature variable
  const ADVariableValue & _temp;

  /// Gradient of temperature variable
  const ADVariableGradient & _grad_temp;

  /// Curvature variable
  const ADVariableValue & _curvature;

  /// Thermal expansion coefficient
  const Real _thermal_expansion;

  /// Reference temperature
  const Real _reference_temperature;

  /// Permeability in Darcy term
  const ADMaterialProperty(Real) & _permeability;

  /// Surface tension coefficient
  const Real & _sigma;

  /// Thermal-capillary coefficient
  const Real & _sigmaT;

  /// Liquid density
  const Real & _rho_l;

  /// Melt pool momentum source
  ADMaterialProperty(RealVectorValue) & _melt_pool_momentum_source;

  /// Laser position
  const MeltPoolLevelSetLocation & _location;

  usingINSMaterialMembers;
};
