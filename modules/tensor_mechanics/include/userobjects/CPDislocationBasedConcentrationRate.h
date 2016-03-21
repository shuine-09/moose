/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CPDISLOCATIONBASEDCONCENTRATIONRATE_H
#define CPDISLOCATIONBASEDCONCENTRATIONRATE_H

#include "CrystalPlasticityStateVarRateComponent.h"
#include "RankTwoTensor.h"

class CPDislocationBasedConcentrationRate;

template<>
InputParameters validParams<CPDislocationBasedConcentrationRate>();

/**
 * A general dislocation based constitutive model userobject class to calculate dislocation density rate component
 * following rho_dot^alpha = k * rho^alpha * abs(gamma_dot^alpha)
 */
class CPDislocationBasedConcentrationRate : public CrystalPlasticityStateVarRateComponent
{
 public:
  CPDislocationBasedConcentrationRate(const InputParameters & parameters);

  virtual bool calcStateVariableEvolutionRateComponent(unsigned int qp, std::vector<Real> & val) const;

 protected:

  const MaterialProperty<std::vector<Real> > & _cv;
  const MaterialProperty<std::vector<RankTwoTensor> > & _flow_direction;
  const MaterialProperty<std::vector<Real> > & _rho_m;
  const MaterialProperty<std::vector<Real> > & _rho_i;
  const MaterialProperty<RankTwoTensor> & _pk2;

  Real _prefactor;

  Real _c0;
  Real _diffusivity;
  Real _rc;
  Real _b;
  Real _molar_volume;
  Real _gas_constant;
  Real _temp;
};

#endif
