//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Kernel.h"
#include "DerivativeMaterialInterface.h"

/**
 * This class computes residual term of damage variable due to pressure.
 */

class PhaseFieldPressurizedFractureDamage : public DerivativeMaterialInterface<Kernel>
{
public:
  static InputParameters validParams();

  PhaseFieldPressurizedFractureDamage(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  /// Number of displacement variables
  unsigned int _ndisp;

  /// Gradient of displacements
  std::vector<const VariableGradient *> _grad_disp;

  /// Gradient of displacements
  std::vector<const VariableValue *> _disp;

  /// Displacement variables IDs
  std::vector<unsigned int> _disp_var;

  /// Material property defining pressure, declared elsewhere
  const MaterialProperty<Real> & _pressure;
  const MaterialProperty<Real> & _L;
};
