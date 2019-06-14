/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "VesicleAreaConstraint.h"

registerMooseObject("MooseApp", VesicleAreaConstraint);

template <>
InputParameters
validParams<VesicleAreaConstraint>()
{
  InputParameters params = validParams<NonlocalKernel>();
  params.addRequiredParam<UserObjectName>("user_object",
                                          "Name of a SimpleTestShapeElementUserObject");
  params.addParam<Real>("alpha_a", 1000, "The penalty parameter of vesicle area.");
  params.addParam<Real>("epsilon", 0.01, "The interfacial penalty parameter.");
  return params;
}

VesicleAreaConstraint::VesicleAreaConstraint(const InputParameters & parameters)
  : NonlocalKernel(parameters),
    _shp(getUserObject<VesicleAreaConstraintUserObject>("user_object")),
    _shp_integral(_shp.getIntegral()),
    _shp_jacobian(_shp.getJacobian()),
    _var_dofs(_var.dofIndices()),
    _alpha_a(getParam<Real>("alpha_a")),
    _epsilon(getParam<Real>("epsilon")),
    _second_u(second())
{
}

Real
VesicleAreaConstraint::computeQpResidual()
{
  Real rz_coord = _q_point[_qp](0);

  return _test[_i][_qp] * _shp_integral * _alpha_a *
         (-3.0 / std::sqrt(2.0) * _epsilon * (_second_u[_qp].tr() + _grad_u[_qp](0) / rz_coord));
}

Real
VesicleAreaConstraint::computeQpJacobian()
{
  Real rz_coord = _q_point[_qp](0);

  return _test[_i][_qp] * _shp_jacobian[_var_dofs[_j]] * _alpha_a *
         (-3.0 / std::sqrt(2.0) * _epsilon * (_second_u[_qp].tr() + _grad_u[_qp](0) / rz_coord));
}

Real
VesicleAreaConstraint::computeQpNonlocalJacobian(dof_id_type dof_index)
{
  Real rz_coord = _q_point[_qp](0);

  return _test[_i][_qp] * _shp_jacobian[dof_index] * _alpha_a *
         (-3.0 / std::sqrt(2.0) * _epsilon * (_second_u[_qp].tr() + _grad_u[_qp](0) / rz_coord));
}
