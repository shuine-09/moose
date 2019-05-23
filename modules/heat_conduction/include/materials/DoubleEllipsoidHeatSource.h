//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"

// Forward Declarations
class DoubleEllipsoidHeatSource;

template <>
InputParameters validParams<DoubleEllipsoidHeatSource>();

/**
 * Double ellipsoid heat source distribution.
 */
class DoubleEllipsoidHeatSource : public Material
{
public:
  DoubleEllipsoidHeatSource(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// power
  const Real _P;
  /// process efficienty
  const Real _eta;
  /// transverse ellipsoid axe
  const Real _a;
  /// depth ellipsoid axe
  const Real _b;
  /// longitudinal ellipsoid axe
  const Real _c;
  /// heating spot travel speed
  const Real _v;
  /// scaling factor
  const Real _f;

  ADMaterialProperty<Real> & _volumetric_heat;
};
