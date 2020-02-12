//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"

class PFCrackOpeningDisplacement;

template <>
InputParameters validParams<PFCrackOpeningDisplacement>();

class PFCrackOpeningDisplacement : public AuxKernel
{
public:
  PFCrackOpeningDisplacement(const InputParameters & parameters);
  virtual ~PFCrackOpeningDisplacement() {}

protected:
  virtual Real computeValue();

  unsigned int _ndisp;
  std::vector<const VariableValue *> _disp;
  const VariableGradient & _grad_c;
};
