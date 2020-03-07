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

template <>
InputParameters validParams<CPDislocationAndPptAthermalResistance>();

class CPDislocationAndPptAthermalResistance : public CPDislocationBasedAthermalSlipResistance
{
public:
  CPDislocationAndPptAthermalResistance(const InputParameters & parameters);

  virtual bool calcSlipResistance(unsigned int qp, std::vector<Real> & val) const;

protected:
  const Function & _precipitate_radius;
  Real _precipitate_volume_fraction;
  Real _orowan_strength_factor;
};

#endif
