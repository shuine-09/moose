/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CPDISLOCATIONBASEDATHERMALSLIPRESISTANCE_H
#define CPDISLOCATIONBASEDATHERMALSLIPRESISTANCE_H

#include "CrystalPlasticitySlipResistance.h"

class CPDislocationBasedAthermalSlipResistance;

template<>
InputParameters validParams<CPDislocationBasedAthermalSlipResistance>();

/**
 * Dislocation based constitutive mode userobject class for athermal slip resistance.
 */
class CPDislocationBasedAthermalSlipResistance : public CrystalPlasticitySlipResistance
{
 public:
  CPDislocationBasedAthermalSlipResistance(const InputParameters & parameters);

  virtual bool calcSlipResistance(unsigned int qp, std::vector<Real> & val) const;

 protected:

  const MaterialProperty<std::vector<Real> > & _mat_prop_mobile_dislocation_density;

  const MaterialProperty<std::vector<Real> > & _mat_prop_immobile_dislocation_density;

  std::vector<Real> _interaction_matrix;

  Real _b;
  Real _q_p;
  Real _shear_mod;
  Real _self_harden;
  Real _latent_harden;
};

#endif // CPDISLOCATIONBASEDATHERMALSLIPRESISTANCE_H
