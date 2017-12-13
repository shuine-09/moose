/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef XFEMRESIDUALOPENING_H
#define XFEMRESIDUALOPENING_H

// MOOSE includes
#include "XFEMMaterialManagerConstraint.h"
#include "MooseMesh.h"
// #include "ComputeCohesiveTraction.h"

// Forward Declarations
class XFEMResidualOpening;

template <>
InputParameters validParams<XFEMResidualOpening>();

class XFEMResidualOpening : public XFEMMaterialManagerConstraint
{
public:
  XFEMResidualOpening(const InputParameters & parameters);
  virtual ~XFEMResidualOpening();

protected:
  /**
   * Set information needed for constraint integration
   */
  virtual void reinitConstraintQuadrature(const ElementPairInfo & element_pair_info) override;
  /**
   *  Compute the residual for one of the constraint quadrature points.
   */
  virtual Real computeQpResidual(Moose::DGResidualType type) override;

  /**
   *  Compute the Jacobian for one of the constraint quadrature points.
   */
  virtual Real computeQpJacobian(Moose::DGJacobianType type) override;

  virtual void initialSetup() override;

  /// Vector normal to the internal interface
  Point _interface_normal;

  const VariableValue & _disp_x;
  const VariableValue & _disp_x_neighbor;
  const VariableValue & _disp_y;
  const VariableValue & _disp_y_neighbor;

  const unsigned int _component;

  const MaterialProperty<Real> * _max_normal_separation_old;
  const MaterialProperty<Real> * _normal_separation;

  const std::string _base_name;

  Real _alpha;

  Real _roughness_fuel;
};

#endif /* XFEMRESIDUALOPENING_H */
