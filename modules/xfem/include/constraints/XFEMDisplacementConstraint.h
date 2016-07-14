/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef XFEMDISPLACEMENTCONSTRAINT_H
#define XFEMDISPLACEMENTCONSTRAINT_H

// MOOSE includes
#include "ElemElemConstraint.h"
#include "MooseMesh.h"

// Forward Declarations
class XFEMDisplacementConstraint;

template<>
InputParameters validParams<XFEMDisplacementConstraint>();

class XFEMDisplacementConstraint : public ElemElemConstraint
{
public:
  XFEMDisplacementConstraint(const InputParameters & parameters);
  virtual ~XFEMDisplacementConstraint();

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

  /// Stabilization parameter in Nitsche's formulation
  Real _alpha;
};

#endif /* XFEMDISPLACEMENTCONSTRAINT_H_ */
