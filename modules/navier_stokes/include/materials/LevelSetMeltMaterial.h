/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#pragma once

#include "ADMaterial.h"

#define usingLevelSetMeltMateriallMembers usingMaterialMembers;

template <ComputeStage>
class LevelSetMeltMaterial;

declareADValidParams(LevelSetMeltMaterial);

/**
 *
 */
template <ComputeStage compute_stage>
class LevelSetMeltMaterial : public ADMaterial<compute_stage>
{
public:
  LevelSetMeltMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// Level set variable
  const ADVariableValue & _ls;

  /// Density
  ADMaterialProperty(Real) & _rho;

  /// Viscosity
  ADMaterialProperty(Real) & _mu;

  /// Gas density
  const Real & _rho_g;

  /// Liquid density
  const Real & _rho_l;

  /// Solid density
  const Real & _rho_s;

  /// Gas viscosity
  const Real & _mu_g;

  /// Liquid viscosity
  const Real & _mu_l;

  /// Solid viscosity
  const Real & _mu_s;

  /// Liquid mass fraction
  const ADMaterialProperty(Real) & _f_l;

  /// Solid mass fraction
  const ADMaterialProperty(Real) & _f_s;

  /// Liquid volume fraction
  const ADMaterialProperty(Real) & _g_l;

  /// Solid volume fraction
  const ADMaterialProperty(Real) & _g_s;

  /// Permeability in Darcy term
  ADMaterialProperty(Real) & _permeability;

  /// Constant in Darcy term
  const Real & _K0;

  usingMaterialMembers;
};
