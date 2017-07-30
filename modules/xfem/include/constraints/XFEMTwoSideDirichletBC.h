/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef XFEMTWOSIDEDIRICHLETBC_H
#define XFEMTWOSIDEDIRICHLETBC_H

// MOOSE includes
#include "ElemElemConstraint.h"
#include "MooseMesh.h"

// Forward Declarations
class XFEMTwoSideDirichletBC;

template <>
InputParameters validParams<XFEMTwoSideDirichletBC>();

class XFEMTwoSideDirichletBC : public ElemElemConstraint
{
public:
  XFEMTwoSideDirichletBC(const InputParameters & parameters);
  virtual ~XFEMTwoSideDirichletBC();

protected:
  /**
   * Set information needed for constraint integration
   */
  virtual void reinitConstraintQuadrature(const ElementPairInfo & element_pair_info);

  /**
   *  Compute the residual for one of the constraint quadrature points.
   */
  virtual Real computeQpResidual(Moose::DGResidualType type);

  /**
   *  Compute the Jacobian for one of the constraint quadrature points.
   */
  virtual Real computeQpJacobian(Moose::DGJacobianType type);

  /// Vector normal to the internal interface
  Point _interface_normal;

  Real _alpha;

  Real _value1;

  Real _value2;
};

#endif /* XFEMTWOSIDEDIRICHLETBC_H */
