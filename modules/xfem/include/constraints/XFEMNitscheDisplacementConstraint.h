/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef XFEMNITSCHEDISPLACEMENTONSTRAINT_H
#define XFEMNITSCHEDISPLACEMENTONSTRAINT_H

// MOOSE includes
#include "XFEMTwoMaterialManagerConstraint.h"
#include "MooseMesh.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

// Forward Declarations
class XFEMNitscheDisplacementConstraint;

template <>
InputParameters validParams<XFEMNitscheDisplacementConstraint>();

class XFEMNitscheDisplacementConstraint : public XFEMTwoMaterialManagerConstraint
{
public:
  XFEMNitscheDisplacementConstraint(const InputParameters & parameters);
  virtual ~XFEMNitscheDisplacementConstraint();

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

  Real _alpha;

  const unsigned int _component;

  const MaterialProperty<RankTwoTensor> * _stress;
  const MaterialProperty<RankTwoTensor> * _stress_neighbor;

  const MaterialProperty<RankFourTensor> * _elasticity_tensor;
  const MaterialProperty<RankFourTensor> * _elasticity_tensor_neighbor;

  const std::string _base_name;
};

#endif /* XFEMNITSCHEDISPLACEMENTONSTRAINT_H */
