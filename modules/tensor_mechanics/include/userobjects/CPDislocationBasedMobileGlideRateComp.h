/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CPDISLOCATIONBASEDMOBILEGLIDERATECOMP_H
#define CPDISLOCATIONBASEDMOBILEGLIDERATECOMP_H

#include "CrystalPlasticityStateVarRateComponent.h"

class CPDislocationBasedMobileGlideRateComp;

template <>
InputParameters validParams<CPDislocationBasedMobileGlideRateComp>();

/**
 * Dislocation based constitutive model userobject class for mobile dislocation density rate
 * component for glide mechanism.
 */
class CPDislocationBasedMobileGlideRateComp : public CrystalPlasticityStateVarRateComponent
{
public:
  CPDislocationBasedMobileGlideRateComp(const InputParameters & parameters);

  virtual bool calcStateVariableEvolutionRateComponent(unsigned int qp,
                                                       std::vector<Real> & val) const;

protected:
  const MaterialProperty<std::vector<Real>> & _mat_prop_mobile_dislocation_density;

  const MaterialProperty<std::vector<Real>> & _mat_prop_immobile_dislocation_density;

  const MaterialProperty<std::vector<Real>> & _mat_prop_glide_slip_rate;

  const MaterialProperty<std::vector<Real>> & _mat_prop_climb_rate;

  Real _b;
  Real _k_mul;
  Real _r_c;
  Real _beta_rho;

  Real _r_precipitate;
  /// radius corresponding to the loss of coherence
  Real _r_cl;

  /// Radius corrsponding the transition from shearable to non-shearable
  Real _r_trans;

  Real _precipitate_radius;

  Real _precipitate_volume_fraction;
};

#endif // CPDISLOCATIONBASEDMOBILEGLIDERATECOMP_H
