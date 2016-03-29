/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CPDISLOCATIONBASEDAPBSLIPRESISTANCE_H
#define CPDISLOCATIONBASEDAPBSLIPRESISTANCE_H

#include "CrystalPlasticitySlipResistance.h"

class CPDislocationBasedAPBSlipResistance;

template<>
InputParameters validParams<CPDislocationBasedAPBSlipResistance>();

/**
 * Dislocation based constitutive mode userobject class for apb shear slip resistance.
 */
class CPDislocationBasedAPBSlipResistance : public CrystalPlasticitySlipResistance
{
 public:
  CPDislocationBasedAPBSlipResistance(const InputParameters & parameters);

  virtual bool calcSlipResistance(unsigned int qp, std::vector<Real> & val) const;

 protected:

  Real _apb_shear_resist; 
};

#endif // CPDISLOCATIONBASEDSLIPRESISTANCE_H
