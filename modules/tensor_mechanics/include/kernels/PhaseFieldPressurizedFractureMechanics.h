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
 * This class computes residual term of displacement variable due to pressure.
 */

class PhaseFieldPressurizedFractureMechanics : public DerivativeMaterialInterface<Kernel>
{
public:
  static InputParameters validParams();

  PhaseFieldPressurizedFractureMechanics(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  /// An integer corresponding to the direction this kernel acts in
  const unsigned int _component;

  /// Coupled order parameter defining the crack
  const VariableGradient & _grad_c;

  const unsigned int _c_var;

  /// Material property defining pressure, declared elsewhere
  const MaterialProperty<Real> & _pressure;
};
