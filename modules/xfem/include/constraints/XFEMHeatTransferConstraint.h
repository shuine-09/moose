/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef XFEMHEATTRANSFERCONSTRAINT_H
#define XFEMHEATTRANSFERCONSTRAINT_H

// MOOSE includes
#include "ElemElemConstraint.h"
#include "MooseMesh.h"
#include "MaterialPropertyInterface.h"

// Forward Declarations
class XFEMHeatTransferConstraint;

template<>
InputParameters validParams<XFEMHeatTransferConstraint>();

class XFEMHeatTransferConstraint : 
  public ElemElemConstraint,
  public MaterialPropertyInterface
{
public:
  XFEMHeatTransferConstraint(const InputParameters & parameters);
  virtual ~XFEMHeatTransferConstraint();

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

  const MaterialProperty<Real> & _heat_flux;
};

#endif /* XFEMHEATTRANSFERCONSTRAINT_H_ */
