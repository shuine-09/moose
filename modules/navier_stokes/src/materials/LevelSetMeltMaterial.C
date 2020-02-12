//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LevelSetMeltMaterial.h"
#include "Function.h"

registerADMooseObject("NavierStokesApp", LevelSetMeltMaterial);

defineADValidParams(LevelSetMeltMaterial,
                    ADMaterial,
                    params.addRequiredCoupledVar("level_set", "Level set variable");
                    params.addRequiredParam<Real>("rho_g", "Gas density.");
                    params.addRequiredParam<Real>("rho_s", "Solid density.");
                    params.addRequiredParam<Real>("rho_l", "Liquid density.");
                    params.addRequiredParam<Real>("mu_g", "Gas viscosity.");
                    params.addRequiredParam<Real>("mu_l", "Liquid viscosity.");
                    params.addRequiredParam<Real>("mu_s", "Solid viscosity.");
                    params.addRequiredParam<Real>("permeability_constant",
                                                  "Permeability constant"););

template <ComputeStage compute_stage>
LevelSetMeltMaterial<compute_stage>::LevelSetMeltMaterial(const InputParameters & parameters)
  : ADMaterial<compute_stage>(parameters),
    _ls(adCoupledValue("level_set")),
    _rho(declareADProperty<Real>("rho")),
    _mu(declareADProperty<Real>("mu")),
    _rho_g(getParam<Real>("rho_g")),
    _rho_l(getParam<Real>("rho_l")),
    _rho_s(getParam<Real>("rho_s")),
    _mu_g(getParam<Real>("mu_g")),
    _mu_l(getParam<Real>("mu_l")),
    _mu_s(getParam<Real>("mu_s")),
    _f_l(getADMaterialProperty<Real>("liquid_mass_fraction")),
    _f_s(getADMaterialProperty<Real>("solid_mass_fraction")),
    _g_l(getADMaterialProperty<Real>("liquid_volume_fraction")),
    _g_s(getADMaterialProperty<Real>("solid_volume_fraction")),
    _permeability(declareADProperty<Real>("permeability")),
    _K0(getParam<Real>("permeability_constant"))
{
}

template <ComputeStage compute_stage>
void
LevelSetMeltMaterial<compute_stage>::computeQpProperties()
{
  ADReal rho_m = _g_s[_qp] * _rho_s + _g_l[_qp] * _rho_l;
  _rho[_qp] = (1 - _ls[_qp]) * rho_m + _ls[_qp] * _rho_g;

  // ADReal mu_m = (_f_s[_qp] * _mu_s + _f_l[_qp] * _mu_l) * (1 - _ls[_qp]);
  ADReal mu_m = _mu_l * rho_m / _rho_l;
  _mu[_qp] = (1 - _ls[_qp]) * mu_m + _ls[_qp] * _mu_g;

  ADReal f_l = _f_l[_qp] * (1 - _ls[_qp]);

  _permeability[_qp] = mu_m / _K0 * Utility::pow<2>(1 - f_l) / (Utility::pow<3>(f_l) + 1.0e-3);
}
