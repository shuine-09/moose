//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef XFEMVelocitySteelOxLimDiff_H
#define XFEMVelocitySteelOxLimDiff_H

#include "XFEMMovingInterfaceVelocityBase.h"

/**
 *
 * Computes the velocity of the "oxide/gas" interface, i.e. oxide (spinel) surface.
 * Length unit has to be nanometer. Time unit has to be hour.
 *
 * Specific to the C4 model for the  corrosion of steel 21-2N at 700C in a CO2
 * atmosphere, with a fixed metal/oxide interface and a single oxide layer
 * growing outward and composed of MnCr2O4 spinel.
 * The variable used is the atomic Mn concentration [at/nm^3].
 */


class XFEMVelocitySteelOxLimDiff;

template <>
InputParameters validParams<XFEMVelocitySteelOxLimDiff>();

class XFEMVelocitySteelOxLimDiff : public XFEMMovingInterfaceVelocityBase
{
public:
  XFEMVelocitySteelOxLimDiff(const InputParameters & parameters);
  virtual ~XFEMVelocitySteelOxLimDiff() {}

  virtual Real computeMovingInterfaceVelocity(unsigned int point_id) const override;


protected:

  /// Diffusivity of Mn in the steel. [nm²/hr]
  Real _diffusivity_Mn_steel;

  // Temperature [K]
  //Real _temperature;
};

#endif
