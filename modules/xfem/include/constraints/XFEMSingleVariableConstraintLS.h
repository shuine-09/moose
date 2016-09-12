/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef XFEMEQUALVALUECONSTRAINTLS_H
#define XFEMEQUALVALUECONSTRAINTLS_H

// MOOSE includes
#include "ElemElemConstraint.h"
#include "MooseMesh.h"

class XFEM;

// Forward Declarations
class XFEMSingleVariableConstraintLS;

template<>
InputParameters validParams<XFEMSingleVariableConstraintLS>();

class XFEMSingleVariableConstraintLS : public ElemElemConstraint
{
public:
  XFEMSingleVariableConstraintLS(const InputParameters & parameters);
  virtual ~XFEMSingleVariableConstraintLS();

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

  virtual bool shouldReverseSign();
  
  /// Vector normal to the internal interface
  Point _interface_normal;

  /// Stabilization parameter in Nitsche's formulation
  Real _alpha;

  /// Vector normal to the internal interface
  Real _jump;

  /// Vector normal to the internal interface
  Real _jump_flux;

  MooseSharedPointer<XFEM> _xfem;

  /// name of level set variable
  NonlinearVariableName _ls_var_name;

  /// level set variable
  MooseVariable & _ls_var;

  SystemBase & _aux_system;
  const NumericVector<Number> * _aux_solution;

};

#endif /* XFEMEQUALVALUECONSTRAINTLS_H_ */
