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

class XFEMC4VelocityZrOxOut;

template <>
InputParameters validParams<XFEMC4VelocityZrOxOut>();

class XFEMC4VelocityZrOxOut : public XFEMMovingInterfaceVelocityBase
{
public:
  XFEMC4VelocityZrOxOut(const InputParameters & parameters);
  virtual ~XFEMC4VelocityZrOxOut() {}

  virtual Real computeMovingInterfaceVelocity(unsigned int point_id) const override;


protected:

  /// Pointer to XFEMC4VelocityZrOxA object
  const XFEMC4VelocityZrOxA * _velocity_uo;

};

#endif // XFEMC4VELOCITYZROXA_H
