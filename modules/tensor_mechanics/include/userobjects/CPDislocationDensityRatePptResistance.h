/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CPDISLOCATIONDENSITYRATEPPTRESISTANCE_H
#define CPDISLOCATIONDENSITYRATEPPTRESISTANCE_H

#include "CrystalPlasticityStateVarRateComponent.h"
#include "Function.h"

class CPDislocationDensityRatePptResistance;

template<>
InputParameters validParams<CPDislocationDensityRatePptResistance>();

class CPDislocationDensityRatePptResistance : public CrystalPlasticityStateVarRateComponent
{
 public:
  CPDislocationDensityRatePptResistance(const InputParameters & parameters);

  virtual bool calcStateVariableEvolutionRateComponent(unsigned int qp, std::vector<Real> & val) const;

 protected:

  const MaterialProperty<std::vector<Real> > & _mat_prop_rate;
  const MaterialProperty<std::vector<Real> > & _mat_prop;

  Real _number_density;
  Real _size;
  Function * _factor_function;
};

#endif
