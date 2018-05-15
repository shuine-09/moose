//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef XFEMTWOSIDEDIRICHLET_H
#define XFEMTWOSIDEDIRICHLET_H

// MOOSE includes
#include "ElemElemConstraint.h"

// Forward Declarations
class XFEMTwoSideDirichlet;

class XFEM;

template <>
InputParameters validParams<XFEMTwoSideDirichlet>();

class XFEMTwoSideDirichlet : public ElemElemConstraint
{
public:
  XFEMTwoSideDirichlet(const InputParameters & parameters);
  virtual ~XFEMTwoSideDirichlet();

protected:
  virtual void reinitConstraintQuadrature(const ElementPairInfo & element_pair_info) override;

  virtual Real computeQpResidual(Moose::DGResidualType type) override;

  virtual Real computeQpJacobian(Moose::DGJacobianType type) override;

  /// Vector normal to the internal interface
  Point _interface_normal;

  // Penalty parameter in penalty's formulation
  Real _alpha;

  /// Value at positive level set interface
  Real _value_at_positive_level_set_interface;

  /// Value at negative level set interface
  Real _value_at_negative_level_set_interface;

  /// Pointer to the XFEM controller object
  std::shared_ptr<XFEM> _xfem;
};

#endif /* XFEMTWOSIDEDIRICHLET_H_ */
