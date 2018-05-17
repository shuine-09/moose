//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BiMaterialTimeDerivative.h"

// MOOSE includes
#include "Assembly.h"
#include "MooseVariableFE.h"

#include "libmesh/quadrature.h"

registerMooseObject("MooseApp", BiMaterialTimeDerivative);

template <>
InputParameters
validParams<BiMaterialTimeDerivative>()
{
  InputParameters params = validParams<TimeDerivative>();
  params.addClassDescription("The time derivative operator with the weak form of $(\\psi_i, "
                             "\\frac{\\partial u_h}{\\partial t})$.");
  return params;
}

BiMaterialTimeDerivative::BiMaterialTimeDerivative(const InputParameters & parameters)
  : TimeDerivative(parameters), _time_step_scale(getMaterialProperty<Real>("time_step_scale"))
{
}

Real
BiMaterialTimeDerivative::computeQpResidual()
{
  return _time_step_scale[_qp] * TimeDerivative::computeQpResidual();
}

Real
BiMaterialTimeDerivative::computeQpJacobian()
{
  return _time_step_scale[_qp] * TimeDerivative::computeQpJacobian();
}
