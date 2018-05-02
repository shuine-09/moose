//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef WELDINGHEATCONDUCTIONMATERIAL_H
#define WELDINGHEATCONDUCTIONMATERIAL_H

#include "Material.h"

// Forward Declarations
class WeldingHeatConductionMaterial;
class WeldStateIndicator;
class Function;

template <>
InputParameters validParams<WeldingHeatConductionMaterial>();

/**
 * Simple material with constant properties.
 */
class WeldingHeatConductionMaterial : public HeatConductionMaterial
{
public:
  WeldingHeatConductionMaterial(const InputParameters & parameters);

protected:
  virtual void computeProperties();

  const WeldStateIndicator * _weld_state;
};

#endif // WELDINGHEATCONDUCTIONMATERIAL_H
