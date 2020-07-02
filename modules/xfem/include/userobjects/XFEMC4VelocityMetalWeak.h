//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef XFEMC4VELOCITYMETALWEAK_H
#define XFEMC4VELOCITYMETALWEAK_H

#include "XFEMMovingInterfaceVelocityBase.h"

class XFEMC4VelocityMetalWeak;

template <>
InputParameters validParams<XFEMC4VelocityMetalWeak>();

class XFEMC4VelocityMetalWeak : public XFEMMovingInterfaceVelocityBase
{
public:
  XFEMC4VelocityMetalWeak(const InputParameters & parameters);
  virtual ~XFEMC4VelocityMetalWeak() {}

  virtual Real computeMovingInterfaceVelocity(unsigned int point_id) const override;

protected:
  /// Diffusivity of oxygen in the Zr alpha phase
  Real _diffusivity_alpha;

  /// Diffusivity of oxygen in the Zr beta phase
  Real _diffusivity_beta;

  /// Initial position of the oxide/metal interface
  Real _x0;
};

#endif // XFEMC4VELOCITYMETALWEAK_H
