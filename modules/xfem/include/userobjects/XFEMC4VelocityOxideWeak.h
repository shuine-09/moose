//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef XFEMC4VELOCITYOXIDEWEAK_H
#define XFEMC4VELOCITYOXIDEWEAK_H

#include "XFEMMovingInterfaceVelocityBase.h"

class XFEMC4VelocityOxideWeak;

template <>
InputParameters validParams<XFEMC4VelocityOxideWeak>();

class XFEMC4VelocityOxideWeak : public XFEMMovingInterfaceVelocityBase
{
public:
  XFEMC4VelocityOxideWeak(const InputParameters & parameters);
  virtual ~XFEMC4VelocityOxideWeak() {}

  virtual Real computeMovingInterfaceVelocity(unsigned int point_id) const override;

protected:

  /// Diffusivity of oxygen in the Zr alpha phase
  Real _diffusivity_alpha;

  // Initial position of the oxide/metal interface
  Real _x0;
};

#endif // XFEMC4VELOCITYOXIDEWEAK_H
