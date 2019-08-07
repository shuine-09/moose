//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef XFEMC4OXIDEVELOCITY_H
#define XFEMC4OXIDEVELOCITY_H

#include "XFEMMovingInterfaceVelocityBase.h"

class XFEMC4OxideVelocity;

template <>
InputParameters validParams<XFEMC4OxideVelocity>();

class XFEMC4OxideVelocity : public XFEMMovingInterfaceVelocityBase
{
public:
  XFEMC4OxideVelocity(const InputParameters & parameters);
  virtual ~XFEMC4OxideVelocity() {}

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

#endif // XFEMC4OXIDEVELOCITY_H
