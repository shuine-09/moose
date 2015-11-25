/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CRYSTALPLASTICITYSTATEVARIABLEEVOLUTIONRATECOMPONENTGSS_H
#define CRYSTALPLASTICITYSTATEVARIABLEEVOLUTIONRATECOMPONENTGSS_H

#include "CrystalPlasticityStateVariableEvolutionRateComponent.h"
#include "RankTwoTensor.h"

class CrystalPlasticityStateVariableEvolutionRateComponentGSS;


template<>InputParameters validParams<CrystalPlasticityStateVariableEvolutionRateComponentGSS>();

/**
 * Crystal plasticity state variable evolution rate component userobject class
 * The virtual functions written below must be
 * over-ridden in derived classes to provide actual values
 */

class CrystalPlasticityStateVariableEvolutionRateComponentGSS : public CrystalPlasticityStateVariableEvolutionRateComponent
{
 public:
  CrystalPlasticityStateVariableEvolutionRateComponentGSS(const InputParameters & parameters);

  virtual bool calcStateVariableEvolutionRateComponent(unsigned int qp, std::vector<Real> & val) const;
 
 protected:

  unsigned int _num_mat_slip_rates;

  std::vector<const MaterialProperty<std::vector<Real> > * > _mat_prop_slip_rates;
};

#endif // CRYSTALPLASTICITYSTATEVARIABLEEVOLUTIONRATECOMPONENTGSS_H
