//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADIntegratedBC.h"

/**
 * Boundary condition for convective heat flux where temperature and heat transfer coefficient are
 * given by material properties.
 */
class ADConvectiveHeatFluxBC : public ADIntegratedBC
{
public:
  static InputParameters validParams();

  ADConvectiveHeatFluxBC(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const Real _heat_convection_coef;
  const Real _ambient_temperature;
  const Real _emissivity;
  const Real _stefan_boltzmann;
  const Real _substrate_thickness;
  const Real _absorptivity;
  const Real _laser_power;
  const Real _laser_radius;
  const Real _laser_velocity;
};
