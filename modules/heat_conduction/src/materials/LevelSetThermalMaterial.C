//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LevelSetThermalMaterial.h"
#include "Function.h"

registerADMooseObject("LevelSetApp", LevelSetThermalMaterial);

defineADValidParams(LevelSetThermalMaterial,
                    ADMaterial,
                    params.addRequiredCoupledVar("level_set", "Level set variable");
                    params.addRequiredCoupledVar("temperature", "Temperature variable");
                    params.addRequiredParam<Real>("c_g", "Gas specific heat.");
                    params.addRequiredParam<Real>("c_s", "Solid specific heat.");
                    params.addRequiredParam<Real>("c_l", "Liquid specific heat.");
                    params.addRequiredParam<Real>("k_g", "Gas heat conductivity.");
                    params.addRequiredParam<Real>("k_s", "Solid conductivity.");
                    params.addRequiredParam<Real>("k_l", "Liquid conductivity.");
                    params.addRequiredParam<Real>("solidus_temperature", "Solidus temperature.");
                    params.addRequiredParam<Real>("latent_heat", "Latent heat."););

template <ComputeStage compute_stage>
LevelSetThermalMaterial<compute_stage>::LevelSetThermalMaterial(const InputParameters & parameters)
  : ADMaterial<compute_stage>(parameters),
    _ls(adCoupledValue("level_set")),
    _temp(adCoupledValue("temperature")),
    _h(declareADProperty<Real>("enthalpy")),
    _k(declareADProperty<Real>("thermal_conductivity")),
    _cp(declareADProperty<Real>("specific_heat")),
    _c_g(getParam<Real>("c_g")),
    _c_s(getParam<Real>("c_s")),
    _c_l(getParam<Real>("c_l")),
    _k_g(getParam<Real>("k_g")),
    _k_s(getParam<Real>("k_s")),
    _k_l(getParam<Real>("k_l")),
    _latent_heat(getParam<Real>("latent_heat")),
    _solidus_temperature(getParam<Real>("solidus_temperature")),
    _f_l(getADMaterialProperty<Real>("liquid_mass_fraction")),
    _f_s(getADMaterialProperty<Real>("solid_mass_fraction")),
    _g_l(getADMaterialProperty<Real>("liquid_volume_fraction")),
    _g_s(getADMaterialProperty<Real>("solid_volume_fraction"))
{
}

template <ComputeStage compute_stage>
void
LevelSetThermalMaterial<compute_stage>::computeQpProperties()
{
  ADReal delta_l = (_c_s - _c_l) * _solidus_temperature + _latent_heat;

  ADReal k_s = 0.0138 * _temp[_qp] + 9.13;

  ADReal f_l = _f_l[_qp] * (1 - _ls[_qp]);
  ADReal f_s = _f_s[_qp] * (1 - _ls[_qp]);
  ADReal g_l = _g_l[_qp];
  ADReal g_s = _g_s[_qp];

  ADReal c_m = (f_s * _c_s + f_l * _c_l) * (1 - _ls[_qp]);
  ADReal k_m = 1.0 / (g_s / k_s + g_l / _k_l);
  ADReal h_m = c_m * _temp[_qp] + f_l * (1 - _ls[_qp]) * delta_l;
  ADReal h_g = _c_g * _temp[_qp];

  _h[_qp] = (1 - _ls[_qp]) * h_m + _ls[_qp] * h_g;
  _k[_qp] = (1 - _ls[_qp]) * k_m + _ls[_qp] * _k_g;

  _cp[_qp] = (1 - _ls[_qp]) * c_m + _ls[_qp] * _c_g;
}
