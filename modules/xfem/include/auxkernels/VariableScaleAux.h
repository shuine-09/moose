//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef VARIABLESCALEAUX_H
#define VARIABLESCALEAUX_H

// MOOSE includes
#include "AuxKernel.h"

// Forward declarations
class VariableScaleAux;

template <>
InputParameters validParams<VariableScaleAux>();

/**
 * Extract a component from the gradient of a variable
 */
class VariableScaleAux : public AuxKernel
{
public:
  /**
   * Class constructor
   * @param parameters Input parameters for the object
   */
  VariableScaleAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

private:
  /// The variable number of the scale variable
  const unsigned int _scale_var_number;

  /// system reference
  const System & _system;

  /// the subproblem solution vector
  const NumericVector<Number> * _solution;

  /// Scale factor
  Real _scale_factor;
};

#endif // VARIABLESCALEAUX_H
