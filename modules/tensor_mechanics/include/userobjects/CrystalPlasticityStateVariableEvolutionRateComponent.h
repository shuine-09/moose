/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CRYSTALPLASTICITYSTATEVARIABLEEVOLUTIONRATECOMPONENT_H
#define CRYSTALPLASTICITYSTATEVARIABLEEVOLUTIONRATECOMPONENT_H

#include "CrystalPlasticityUOBase.h"
#include "RankTwoTensor.h"

class CrystalPlasticityStateVariableEvolutionRateComponent;


template<>InputParameters validParams<CrystalPlasticityStateVariableEvolutionRateComponent>();

/**
 * Crystal plasticity state variable evolution rate component userobject base class
 * The virtual functions written below must be
 * over-ridden in derived classes to provide actual values
 */

class CrystalPlasticityStateVariableEvolutionRateComponent : public CrystalPlasticityUOBase
{
 public:
  CrystalPlasticityStateVariableEvolutionRateComponent(const InputParameters & parameters);

  /// Returns the slip rate  name  (eg "DislocationGlide")
  virtual std::string crystalPlasticityUOName() const;

  virtual bool calcStateVariableEvolutionRateComponent(unsigned int qp, std::vector<Real> & val) const = 0;
 
 protected:
  unsigned int _nss;
};

#endif // CRYSTALPLASTICITYSTATEVARIABLEEVOLUTIONRATECOMPONENT_H
