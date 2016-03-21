/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CPDISLOCATIONANDPPTATHERMALRESISTANCE_H
#define CPDISLOCATIONANDPPTATHERMALRESISTANCE_H

#include "CPDislocationBasedAthermalSlipResistance.h"
#include "Function.h"

class CPDislocationAndPptAthermalResistance;

template<>
InputParameters validParams<CPDislocationAndPptAthermalResistance>();

class CPDislocationAndPptAthermalResistance : public CPDislocationBasedAthermalSlipResistance
{
public:
  CPDislocationAndPptAthermalResistance(const InputParameters & parameters);

  virtual bool calcSlipResistance(unsigned int qp, std::vector<Real> & val) const;

protected:

  const MaterialProperty<std::vector<Real> > & _shear_rate_prop;

  Real _apb_shear_energy;
  Real _number_density;
  Real _size;
  Function * _factor_function;
};

#endif
