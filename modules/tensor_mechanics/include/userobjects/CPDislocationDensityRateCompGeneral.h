/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CPDISLOCATIONDENSITYRATECOMPGENERAL_H
#define CPDISLOCATIONDENSITYRATECOMPGENERAL_H

#include "CrystalPlasticityStateVarRateComponent.h"

class CPDislocationDensityRateCompGeneral;

template<>
InputParameters validParams<CPDislocationDensityRateCompGeneral>();

/**
 * A general dislocation based constitutive model userobject class to calculate dislocation density rate component
 * following rho_dot^alpha = k * rho^alpha * abs(gamma_dot^alpha)
 */
class CPDislocationDensityRateCompGeneral : public CrystalPlasticityStateVarRateComponent
{
 public:
  CPDislocationDensityRateCompGeneral(const InputParameters & parameters);

  virtual bool calcStateVariableEvolutionRateComponent(unsigned int qp, std::vector<Real> & val) const;

 protected:

  const MaterialProperty<std::vector<Real> > & _mat_prop_dislocation_density;
 
  const MaterialProperty<std::vector<Real> > & _mat_prop_rate;

  Real _prefactor;
};

#endif // CPDISLOCATIONDENSITYRATECOMPGENERAL_H
