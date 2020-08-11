//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef XFEMC4VELOCITYOXIDESPINEL_H
#define XFEMC4VELOCITYOXIDESPINEL_H

#include "XFEMMovingInterfaceVelocityBase.h"

class XFEMC4VelocityOxideSpinel;

template <>
InputParameters validParams<XFEMC4VelocityOxideSpinel>();

class XFEMC4VelocityOxideSpinel : public XFEMMovingInterfaceVelocityBase
{
public:
  XFEMC4VelocityOxideSpinel(const InputParameters & parameters);
  virtual ~XFEMC4VelocityOxideSpinel() {}

  virtual Real computeMovingInterfaceVelocity(unsigned int point_id) const override;


protected:

  /// Diffusivity of oxygen in the Zr alpha phase
  //Real _diffusivity_alpha;

  // Temperature [K]
  //Real _temperature;
};

#endif
