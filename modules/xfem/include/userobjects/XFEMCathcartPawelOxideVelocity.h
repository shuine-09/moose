//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef XFEMCATHCARTPAWELOXIDEVELOCITY_H
#define XFEMCATHCARTPAWELOXIDEVELOCITY_H

#include "XFEMMovingInterfaceVelocityBase.h"

class XFEMCathcartPawelOxideVelocity;

template <>
InputParameters validParams<XFEMCathcartPawelOxideVelocity>();

class XFEMCathcartPawelOxideVelocity : public XFEMMovingInterfaceVelocityBase
{
public:
  XFEMCathcartPawelOxideVelocity(const InputParameters & parameters);
  virtual ~XFEMCathcartPawelOxideVelocity() {}

  virtual Real computeMovingInterfaceVelocity(unsigned int point_id) const override;

protected:
  /// Diffusivity in the positive level set region
  Real _diffusivity_at_positive_level_set;

  /// Diffusivity in the negative level set region
  Real _diffusivity_at_negative_level_set;

  /// Jump of the equilibrium concentrations at phase boundary
  Real _equilibrium_concentration_jump;
};

#endif // XFEMCATHCARTPAWELOXIDEVELOCITY_H
