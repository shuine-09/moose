/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "ConcentrationDiffusion.h"
#include "MooseMesh.h"

registerMooseObject("XFEMApp", ConcentrationDiffusion);

template <>
InputParameters
validParams<ConcentrationDiffusion>()
{
  InputParameters params = validParams<Diffusion>();
  params.addClassDescription(
      "Computes residual/Jacobian contribution for diffusion equation with diffusivity.");
  params.addParam<std::string>("base_name",
                               "Optional parameter that allows the user to define "
                               "mutliple materials systems on the same block.");
  return params;
}

ConcentrationDiffusion::ConcentrationDiffusion(const InputParameters & parameters)
  : Diffusion(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _diffusion_coefficient(getMaterialProperty<Real>(_base_name + "diffusion_coefficient")),
    _diffusion_coefficient_dT(hasMaterialProperty<Real>(_base_name + "diffusion_coefficient_dT")
                                  ? &getMaterialProperty<Real>("diffusion_coefficient_dT")
                                  : NULL)
{
}

Real
ConcentrationDiffusion::computeQpResidual()
{
  return _diffusion_coefficient[_qp] * Diffusion::computeQpResidual();
}

Real
ConcentrationDiffusion::computeQpJacobian()
{
  Real jac = _diffusion_coefficient[_qp] * Diffusion::computeQpJacobian();
  if (_diffusion_coefficient_dT)
    jac += (*_diffusion_coefficient_dT)[_qp] * _phi[_j][_qp] * Diffusion::computeQpResidual();
  return jac;
}
