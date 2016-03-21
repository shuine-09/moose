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

template<>InputParameters validParams<CPDislocationBasedMobileGlideRateComp>();

/**
 * Dislocation based constitutive model userobject class for mobile dislocation density rate component for glide mechanism.
 */
class CPDislocationBasedMobileGlideRateComp : public CrystalPlasticityStateVarRateComponent
{
 public:
  CPDislocationBasedMobileGlideRateComp(const InputParameters & parameters);

  virtual bool calcStateVariableEvolutionRateComponent(unsigned int qp, std::vector<Real> & val) const;

 protected:

  const MaterialProperty<std::vector<Real> > & _mat_prop_mobile_dislocation_density;

  const MaterialProperty<std::vector<Real> > & _mat_prop_immobile_dislocation_density;
 
  const MaterialProperty<std::vector<Real> > & _mat_prop_glide_slip_rate;

  Real _b;
  Real _k_mul;
  Real _r_c;
  Real _beta_rho;
};

#endif // CPDISLOCATIONBASEDMOBILEGLIDERATECOMP_H
