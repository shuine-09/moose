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

// Forward Declarations
class WeakToStrongDiscontinuityC4ZrAux;

/**
 * Computes the value of the original variable of the strong discontinuity
 * problem from the weak discontinuity variable.
 * !!! Not working correctly. Only adds the discontinuity value once an element
 * has been fully crossed by the interface. !!!
 */
class WeakToStrongDiscontinuityC4ZrAux : public AuxKernel
{
public:
  static InputParameters validParams();

  WeakToStrongDiscontinuityC4ZrAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  const Real _jump_value;

  /// The variable number of the weak discontinuity variable we solve for
  const unsigned int _weak_variable_number;

  /// The variable number of the level set variable we are operating on
  const unsigned int _level_set_var_number;

  /// shared pointer to XFEM
  std::shared_ptr<XFEM> _xfem;

  /// System reference
  const System & _system;

  const System & _system2;

  /// The subproblem solution vector
  const NumericVector<Number> * _solution;
  /// The subproblem solution vector
  const NumericVector<Number> * _solution2;

  /// add the discontinuity to the weak problem concentration if in the oxide (positive level set)
  bool _add_jump;


};
