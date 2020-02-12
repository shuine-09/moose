/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#pragma once

#include "ADMaterial.h"
#include "DEDLevelSetLocation.h"

#define usingDEDMateriallMembers usingMaterialMembers;

template <ComputeStage>
class DEDMaterial;

declareADValidParams(DEDMaterial);

template <ComputeStage compute_stage>
class DEDMaterial : public ADMaterial<compute_stage>
{
public:
  DEDMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  const ADVectorVariableValue & _velocity;
  const ADVariableValue & _ls;
  const ADVariableValue & _temp;
  const ADVariableGradient & _grad_temp;
  const ADVariableGradient & _grad_ls1;
  const ADVectorVariableValue & _grad_ls;
  const ADVariableValue & _curvature;

  ADMaterialProperty(Real) & _rho;
  ADMaterialProperty(RealVectorValue) & _grad_rho;
  ADMaterialProperty(Real) & _mu;
  ADMaterialProperty(Real) & _h;
  ADMaterialProperty(Real) & _k;
  ADMaterialProperty(Real) & _dhdT;

  ADMaterialProperty(Real) & _rho_m;
  ADMaterialProperty(Real) & _h_m;
  ADMaterialProperty(Real) & _mu_m;
  ADMaterialProperty(Real) & _powder_feed;
  ADMaterialProperty(Real) & _mass_change;
  ADMaterialProperty(RealVectorValue) & _ded_momentum;
  ADMaterialProperty(RealVectorValue) & _grads_T;
  ADMaterialProperty(Real) & _permeability;

  RealVectorValue _gravity;
  Real _thermal_expansion;
  Real _reference_temperature;

  ADMaterialProperty(Real) & _f_l;
  ADMaterialProperty(Real) & _f_s;

  const Real & _rho_g;
  const Real & _rho_s;
  const Real & _rho_l;

  const Real & _mu_g;
  const Real & _mu_l;
  const Real & _mu_s;

  const Real & _c_g;
  const Real & _c_s;
  const Real & _c_l;

  const Real & _k_g;
  const Real & _k_s;
  const Real & _k_l;

  const Real & _solidus_temperature;
  const Real & _liquidus_temperature;
  const Real & _latent_heat;

  const Real & _K0;
  const Real _sigma;
  const Real _sigmaT;

  // const Function & _laser_center_x;
  // const Function & _laser_center_y;
  // const Function & _laser_center_z;

  const DEDLevelSetLocation & _location;

  const Real & _power;
  const Real & _alpha;
  const Real & _Rb;
  const Real & _Ah;
  const Real & _stefan_boltzmann;
  const Real & _varepsilon;
  const Real & _T0;
  ADMaterialProperty(Real) & _heat_source;

  const Real & _mass_rate;
  const Real & _mass_radius;

  usingMaterialMembers;
};
