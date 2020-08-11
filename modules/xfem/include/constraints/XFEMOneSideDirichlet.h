//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef XFEMONESIDEDIRICHLET_H
#define XFEMONESIDEDIRICHLET_H

// MOOSE includes
#include "ElemElemConstraint.h"

// Forward Declarations
class XFEMOneSideDirichlet;

class XFEM;

template <>
InputParameters validParams<XFEMOneSideDirichlet>();

class XFEMOneSideDirichlet : public ElemElemConstraint
{
public:
  XFEMOneSideDirichlet(const InputParameters & parameters);
  virtual ~XFEMOneSideDirichlet();

protected:
  virtual void reinitConstraintQuadrature(const ElementPairInfo & element_pair_info) override;

  virtual Real computeQpResidual(Moose::DGResidualType type) override;

  virtual Real computeQpJacobian(Moose::DGJacobianType type) override;

  /// Vector normal to the internal interface
  Point _interface_normal;

  // Penalty parameter in penalty's formulation
  Real _alpha;

  /// Value of the variable on the constrained side
  Real _value;

  // Boolean specifying the side of the interface the constraint is applied to (positive if true)
  bool _positive_side;

  /// Pointer to the XFEM controller object
  std::shared_ptr<XFEM> _xfem;
};

#endif /* XFEMONESIDEDIRICHLET_H_ */
