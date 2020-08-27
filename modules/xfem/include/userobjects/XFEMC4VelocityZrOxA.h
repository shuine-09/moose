//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef XFEMC4VelocityZrOxA_H
#define XFEMC4VelocityZrOxA_H

#include "XFEMMovingInterfaceVelocityBase.h"


/**
 *
 * Computes the velocity of the oxide/alpha interface.
 * Unit length has to be micrometer.
 *
 * Specific to the weak discontinuity equivalent of the C4 model for the
 * high-temperature corrosion of Zircaloy-4 (1000C to 1500C)
 * The variable used is the weak (atomic) oxygen concentration scaled by
 * the atomic concentration of Zr in the metal.
 */

class XFEMC4VelocityZrOxA;

template <>
InputParameters validParams<XFEMC4VelocityZrOxA>();

class XFEMC4VelocityZrOxA : public XFEMMovingInterfaceVelocityBase
{
public:
  XFEMC4VelocityZrOxA(const InputParameters & parameters);
  virtual ~XFEMC4VelocityZrOxA() {}

  virtual Real computeMovingInterfaceVelocity(unsigned int point_id) const override;


protected:

  /// Diffusivity of oxygen in the Zr alpha phase
  //Real _diffusivity_alpha;

  // Temperature [K]
  Real _temperature;
};

#endif // XFEMC4VELOCITYZROXA_H
