/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CPDISLOCATIONBASEDIMMOBILEDENSITY_H
#define CPDISLOCATIONBASEDIMMOBILEDENSITY_H

#include "CrystalPlasticityStateVariable.h"
#include "CrystalPlasticitySlipResistance.h"

class CPDislocationBasedImmobileDensity;

template<>
InputParameters validParams<CPDislocationBasedImmobileDensity>();

/**
 * Dislocation based constitutive mode userobject class for immobile dislocation density.
 */
class CPDislocationBasedImmobileDensity : public CrystalPlasticityStateVariable
{
 public:
   CPDislocationBasedImmobileDensity(const InputParameters & parameters);

   virtual void initSlipSysProps(std::vector<Real> & val) const;

 protected:
   virtual void assignSlipSysRes(std::vector<Real> & val) const;

   /// Read initial slip system resistances  from .i file
   virtual void getInitSlipSysRes(std::vector<Real> & val) const;

   Real _initial_value;
};

#endif // CPDISLOCATIONBASEDIMMOBILEDENSITY_H
