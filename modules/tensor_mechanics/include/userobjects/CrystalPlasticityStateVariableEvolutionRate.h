/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CRYSTALPLASTICITYSTATEVARIABLEEVOLUTIONRATE_H
#define CRYSTALPLASTICITYSTATEVARIABLEEVOLUTIONRATE_H

#include "CrystalPlasticityUOBase.h"
#include "RankTwoTensor.h"

class CrystalPlasticityStateVariableEvolutionRate;

template<>
InputParameters validParams<CrystalPlasticityStateVariableEvolutionRate>();

/**
 * Crystal plasticity slip rate userobject class
 * The virtual functions written below must be
 * over-ridden in derived classes to provide actual values
 */

class CrystalPlasticityStateVariableEvolutionRate : public CrystalPlasticityUOBase
{
 public:
  CrystalPlasticityStateVariableEvolutionRate(const InputParameters & parameters);

  /// Returns the slip rate  name  (eg "DislocationGlide")
  virtual std::string crystalPlasticityUOName() const;

  virtual bool calcStateVariableEvolutionRate(unsigned int qp, std::vector<Real> & val) const;
 
 protected:

  unsigned int _nss;

  std::vector<const MaterialProperty<std::vector<Real> > * > _mat_prop_state_var_evol_rate_comps;

  unsigned int _num_mat_state_var_evol_rate_comps;

};

#endif // CRYSTALPLASTICITYSTATEVARIABLEEVOLUTIONRATE_H
