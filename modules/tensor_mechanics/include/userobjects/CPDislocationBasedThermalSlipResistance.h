/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CPDISLOCATIONBASEDTHERMALSLIPRESISTANCE_H
#define CPDISLOCATIONBASEDTHERMALSLIPRESISTANCE_H

#include "CrystalPlasticitySlipResistance.h"

class CPDislocationBasedThermalSlipResistance;

template<>
InputParameters validParams<CPDislocationBasedThermalSlipResistance>();

/**
 * Dislocation based constitutive mode userobject class for thermal slip resistance.
 */
class CPDislocationBasedThermalSlipResistance : public CrystalPlasticitySlipResistance
{
 public:
  CPDislocationBasedThermalSlipResistance(const InputParameters & parameters);

  virtual bool calcSlipResistance(unsigned int qp, std::vector<Real> & val) const;

 protected:

  Real _thermal_resist; 
};

#endif // CPDISLOCATIONBASEDTHERMALSLIPRESISTANCE_H
