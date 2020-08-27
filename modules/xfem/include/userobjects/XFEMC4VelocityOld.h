//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef XFEMC4VelocityOld_H
#define XFEMC4VelocityOld_H

#include "XFEMMovingInterfaceVelocityBase.h"

/**
*
* Computes the oxide/alpha interface velocity for Zr C4 model.
* !!! Old Leo's version !!!
*
*/



class XFEMC4VelocityOld;

template <>
InputParameters validParams<XFEMC4VelocityOld>();

class XFEMC4VelocityOld : public XFEMMovingInterfaceVelocityBase
{
public:
  XFEMC4VelocityOld(const InputParameters & parameters);
  virtual ~XFEMC4VelocityOld() {}

  virtual Real computeMovingInterfaceVelocity(unsigned int point_id) const override;

protected:
  /// Diffusivity in the positive level set region
  Real _diffusivity_at_positive_level_set;

  /// Diffusivity in the negative level set region
  Real _diffusivity_at_negative_level_set;

  /// Jump of the equilibrium concentrations at phase boundary
  Real _equilibrium_concentration_jump;

  Real _x0;
};

#endif // XFEMC4VelocityOld_H
