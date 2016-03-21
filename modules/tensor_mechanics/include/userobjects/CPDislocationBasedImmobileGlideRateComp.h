/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CPDISLOCATIONBASEDIMMOBILEGLIDERATECOMP_H
#define CPDISLOCATIONBASEDIMMOBILEGLIDERATECOMP_H

#include "CrystalPlasticityStateVarRateComponent.h"

class CPDislocationBasedImmobileGlideRateComp;

template<>InputParameters validParams<CPDislocationBasedImmobileGlideRateComp>();

/**
 * Dislocation based constitutive model userobject class for immobile dislocation density rate component for glide mechanism.
 */
class CPDislocationBasedImmobileGlideRateComp : public CrystalPlasticityStateVarRateComponent
{
 public:
  CPDislocationBasedImmobileGlideRateComp(const InputParameters & parameters);

  virtual bool calcStateVariableEvolutionRateComponent(unsigned int qp, std::vector<Real> & val) const;

 protected:

  const MaterialProperty<std::vector<Real> > & _mat_prop_mobile_dislocation_density;

  const MaterialProperty<std::vector<Real> > & _mat_prop_immobile_dislocation_density;

  const MaterialProperty<std::vector<Real> > & _mat_prop_glide_slip_rate;

  Real _b;
  Real _beta_rho;
};

#endif // CPDISLOCATIONBASEDIMMOBILEGLIDERATECOMP_H
