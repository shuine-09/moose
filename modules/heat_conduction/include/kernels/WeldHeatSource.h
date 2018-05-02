//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef WELDHEATSOURCE_H
#define WELDHEATSOURCE_H

#include "BodyForce.h"

// Forward Declarations
class WeldHeatSource;
class WeldStateIndicator;

template <>
InputParameters validParams<WeldHeatSource>();

class WeldHeatSource : public BodyForce
{
public:
  WeldHeatSource(const InputParameters & parameters);

  virtual Real computeQpResidual() override;

private:
  ///
  const WeldStateIndicator * _weld_state;
};

#endif
