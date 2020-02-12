//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MushyZoneMaterial.h"
#include "Function.h"

registerADMooseObject("LevelSetApp", MushyZoneMaterial);

defineADValidParams(MushyZoneMaterial,
                    ADMaterial,
                    params.addRequiredCoupledVar("temperature", "Temperature variable");
                    params.addRequiredParam<Real>("rho_s", "Solid density.");
                    params.addRequiredParam<Real>("rho_l", "Liquid density.");
                    params.addRequiredParam<Real>("solidus_temperature", "Solidus temperature.");
                    params.addRequiredParam<Real>("liquidus_temperature",
                                                  "Liquidus temperature."););

template <ComputeStage compute_stage>
MushyZoneMaterial<compute_stage>::MushyZoneMaterial(const InputParameters & parameters)
  : ADMaterial<compute_stage>(parameters),
    _temp(adCoupledValue("temperature")),
    _solidus_temperature(getParam<Real>("solidus_temperature")),
    _liquidus_temperature(getParam<Real>("liquidus_temperature")),
    _f_l(declareADProperty<Real>("liquid_mass_fraction")),
    _f_s(declareADProperty<Real>("solid_mass_fraction")),
    _g_l(declareADProperty<Real>("liquid_volume_fraction")),
    _g_s(declareADProperty<Real>("solid_volume_fraction")),
    _rho_s(getParam<Real>("rho_s")),
    _rho_l(getParam<Real>("rho_l"))
{
}

template <ComputeStage compute_stage>
void
MushyZoneMaterial<compute_stage>::computeQpProperties()
{
  _f_l[_qp] = 1;

  if (_temp[_qp] < _solidus_temperature)
    _f_l[_qp] = 0;
  else if (_temp[_qp] >= _solidus_temperature && _temp[_qp] <= _liquidus_temperature)
    _f_l[_qp] =
        (_temp[_qp] - _solidus_temperature) / (_liquidus_temperature - _solidus_temperature);

  _f_s[_qp] = 1.0 - _f_l[_qp];

  _g_l[_qp] = _f_l[_qp] / _rho_l / ((1 - _f_l[_qp]) / _rho_s + _f_l[_qp] / _rho_l);
  _g_s[_qp] = 1.0 - _g_l[_qp];
}
