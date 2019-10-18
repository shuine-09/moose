//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HeatConduction.h"
#include "MooseMesh.h"

registerMooseObjectAliased("HeatConductionApp", HeatConductionKernel, "HeatConduction");

template <>
InputParameters
validParams<HeatConductionKernel>()
{
  InputParameters params = validParams<Diffusion>();
  params.addClassDescription(
      "Computes residual/Jacobian contribution for $(k \\nabla T, \\nabla \\psi)$ term.");
  params.addParam<MaterialPropertyName>(
      "diffusion_coefficient",
      "thermal_conductivity",
      "Property name of the diffusivity (Default: thermal_conductivity)");
  params.addParam<MaterialPropertyName>(
      "diffusion_coefficient_dT",
      "thermal_conductivity_dT",
      "Property name of the derivative of the diffusivity with respect "
      "to the variable (Default: thermal_conductivity_dT)");
  params.addParam<bool>("coupled_to_damage", false, "Coupled to damage.");
  params.addCoupledVar("c", "Phase field damage variable");
  params.addParam<Real>("kdamage", 0.0, "Stiffness of damaged matrix");
  params.set<bool>("use_displaced_mesh") = true;
  return params;
}

HeatConductionKernel::HeatConductionKernel(const InputParameters & parameters)
  : Diffusion(parameters),
    _diffusion_coefficient(getMaterialProperty<Real>("diffusion_coefficient")),
    _diffusion_coefficient_dT(hasMaterialProperty<Real>("diffusion_coefficient_dT")
                                  ? &getMaterialProperty<Real>("diffusion_coefficient_dT")
                                  : NULL),
    _coupled_to_damage(getParam<bool>("coupled_to_damage")),
    _c(coupledValue("c")),
    _c_var(coupled("c")),
    _kdamage(getParam<Real>("kdamage"))
{
}

Real
HeatConductionKernel::computeQpResidual()
{
  if (_coupled_to_damage)
    return _diffusion_coefficient[_qp] *
           ((1.0 - _c[_qp]) * (1.0 - _c[_qp]) * (1 - _kdamage) + _kdamage) *
           Diffusion::computeQpResidual();
  else
    return _diffusion_coefficient[_qp] * Diffusion::computeQpResidual();
}

Real
HeatConductionKernel::computeQpJacobian()
{
  Real jac = _diffusion_coefficient[_qp] * Diffusion::computeQpJacobian();

  if (_coupled_to_damage)
    jac = _diffusion_coefficient[_qp] * Diffusion::computeQpJacobian() *
          ((1.0 - _c[_qp]) * (1.0 - _c[_qp]) * (1 - _kdamage) + _kdamage);

  if (_diffusion_coefficient_dT)
    jac += (*_diffusion_coefficient_dT)[_qp] * _phi[_j][_qp] * Diffusion::computeQpResidual();
  return jac;
}

Real
HeatConductionKernel::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_coupled_to_damage && jvar == _c_var)
  {
    return _diffusion_coefficient[_qp] * Diffusion::computeQpResidual() * (-2.0 * _phi[_j][_qp]) *
           (1.0 - _c[_qp]) * (1 - _kdamage);
  }

  // Returns if coupled variable is not c (damage variable)
  return 0.0;
}
