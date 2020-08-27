//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ElemElemConstraint.h"

// Forward Declarations
class XFEM;

/**
 * Sets variable value at the oxide/alpha interface.
 *
 * Specific to the weak discontinuity equivalent of the C4 model for the
 * high-temperature corrosion of Zircaloy-4 (1000C to 1500C)
 * The variable used is the weak (atomic) oxygen concentration scaled by
 * the atomic concentration of Zr in the metal.
 */

class XFEMEqualValueAtInterfaceC4aox : public ElemElemConstraint
{
public:
  static InputParameters validParams();

  XFEMEqualValueAtInterfaceC4aox(const InputParameters & parameters);
  virtual ~XFEMEqualValueAtInterfaceC4aox();

protected:
  virtual void reinitConstraintQuadrature(const ElementPairInfo & element_pair_info) override;

  virtual Real computeQpResidual(Moose::DGResidualType type) override;

  virtual Real computeQpJacobian(Moose::DGJacobianType type) override;

  // Penalty parameter in penalty formulation
  Real _alpha;

  /// Value at the interface
  Real _value;

  /// Temperature [K] that will determine the value at the alpha/beta interface
  Real _temperature;

  /// Boolean that specifies whether or not an offset is needed due to transform to weak discontinuity,
  /// i.e. whether or not there is another interface (alpha beta) before this one
  bool _offset;

  /// Pointer to the XFEM controller object
  std::shared_ptr<XFEM> _xfem;
};
