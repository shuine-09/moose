/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CPDISLOCATIONBASEDDENSITYOROWANLOOPINGRATECOMP_H
#define CPDISLOCATIONBASEDDENSITYOROWANLOOPINGRATECOMP_H

#include "CrystalPlasticityStateVarRateComponent.h"

class CPDislocationBasedDensityOrowanLoopingRateComp;

template <>
InputParameters validParams<CPDislocationBasedDensityOrowanLoopingRateComp>();

/**
 * Dislocation based constitutive model userobject class for dislocation density rate
 * component for orowan looping mechanism.
 */
class CPDislocationBasedDensityOrowanLoopingRateComp : public CrystalPlasticityStateVarRateComponent
{
public:
  CPDislocationBasedDensityOrowanLoopingRateComp(const InputParameters & parameters);

  virtual bool calcStateVariableEvolutionRateComponent(unsigned int qp,
                                                       std::vector<Real> & val) const;

protected:
  const MaterialProperty<std::vector<Real>> & _mat_prop_glide_slip_rate;

  /// radius corresponding to the loss of coherence
  Real _r_cl;

  /// Radius corrsponding the transition from shearable to non-shearable
  Real _r_trans;

  Real _factor;

  Real _precipitate_radius;

  Real _precipitate_volume_fraction;
};

#endif // CPDISLOCATIONBASEDDENSITYOROWANLOOPINGRATECOMP_H
