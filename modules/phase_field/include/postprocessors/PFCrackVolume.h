//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ElementIntegralPostprocessor.h"
#include "MooseVariableInterface.h"

// Forward Declarations

/**
 * Compute a volume integral of the phase-field cracks.
 */
class PFCrackVolume : public ElementIntegralPostprocessor, public MooseVariableInterface<Real>
{
public:
  static InputParameters validParams();

  PFCrackVolume(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral();

  MooseVariable & _var;

  /// Holds the solution at current quadrature points
  const VariableValue & _u;

  /// Holds the solution gradient at the current quadrature points
  const VariableGradient & _grad_u;

  /// Holds the solution derivative at the current quadrature points
  const VariableValue & _u_dot;

  /// l
  const Real _l;
};
