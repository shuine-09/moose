/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CPDISLOCATIONBASEDMOBILEDENSITY_H
#define CPDISLOCATIONBASEDMOBILEDENSITY_H

#include "CrystalPlasticityStateVariable.h"
#include "CrystalPlasticitySlipResistance.h"

class CPDislocationBasedMobileDensity;

template<>
InputParameters validParams<CPDislocationBasedMobileDensity>();

/**
 * Dislocation based constitutive mode userobject class for mobile dislocation density.
 */
class CPDislocationBasedMobileDensity : public CrystalPlasticityStateVariable
{
 public:
   CPDislocationBasedMobileDensity(const InputParameters & parameters);

   virtual void initSlipSysProps(std::vector<Real> & val) const;

 protected:
   virtual void assignSlipSysRes(std::vector<Real> & val) const;

   /// Read initial slip system resistances  from .i file
   virtual void getInitSlipSysRes(std::vector<Real> & val) const;

   Real _initial_value;
};

#endif // CPDISLOCATIONBASEDMOBILEDENSITY_H
