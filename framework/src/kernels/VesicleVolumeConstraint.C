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

#include "VesicleVolumeConstraint.h"

registerMooseObject("MooseApp", VesicleVolumeConstraint);

template <>
InputParameters
validParams<VesicleVolumeConstraint>()
{
  InputParameters params = validParams<NonlocalKernel>();
  params.addRequiredParam<UserObjectName>("user_object",
                                          "Name of a SimpleTestShapeElementUserObject");
  params.addParam<Real>("alpha_v", 1000, "The penalty parameter of vesicle volume.");
  return params;
}

VesicleVolumeConstraint::VesicleVolumeConstraint(const InputParameters & parameters)
  : NonlocalKernel(parameters),
    _shp(getUserObject<VesicleVolumeConstraintUserObject>("user_object")),
    _shp_integral(_shp.getIntegral()),
    _shp_jacobian(_shp.getJacobian()),
    _var_dofs(_var.dofIndices()),
    _alpha_v(getParam<Real>("alpha_v"))
{
}

Real
VesicleVolumeConstraint::computeQpResidual()
{
  return _test[_i][_qp] * _shp_integral * _alpha_v;
}

Real
VesicleVolumeConstraint::computeQpJacobian()
{
  return _test[_i][_qp] * _shp_jacobian[_var_dofs[_j]] * _alpha_v;
}

Real
VesicleVolumeConstraint::computeQpNonlocalJacobian(dof_id_type dof_index)
{
  return _test[_i][_qp] * _shp_jacobian[dof_index] * _alpha_v;
}
