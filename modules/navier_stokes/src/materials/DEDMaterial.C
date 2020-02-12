/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "DEDMaterial.h"
#include "Function.h"

registerADMooseObject("NavierStokesApp", DEDMaterial);

defineADValidParams(
    DEDMaterial, ADMaterial, params.addRequiredCoupledVar("level_set", "Level set variable");
    params.addRequiredCoupledVar("velocity", "The velocity");
    params.addRequiredCoupledVar("grad_level_set", "Regularized gradient of Level set variable");
    params.addRequiredCoupledVar("temperature", "Temperature variable");
    params.addRequiredCoupledVar("curvature", "Level set variable");
    params.addRequiredParam<Real>("rho_g", "Gas density.");
    params.addRequiredParam<Real>("rho_s", "Solid density.");
    params.addRequiredParam<Real>("rho_l", "Liquid density.");
    params.addRequiredParam<Real>("mu_g", "Gas viscosity.");
    params.addRequiredParam<Real>("mu_l", "Liquid viscosity.");
    params.addRequiredParam<Real>("mu_s", "Solid viscosity.");
    params.addRequiredParam<Real>("c_g", "Gas specific heat.");
    params.addRequiredParam<Real>("c_s", "Solid specific heat.");
    params.addRequiredParam<Real>("c_l", "Liquid specific heat.");
    params.addRequiredParam<Real>("k_g", "Gas heat conductivity.");
    params.addRequiredParam<Real>("k_s", "Solid conductivity.");
    params.addRequiredParam<Real>("k_l", "Liquid conductivity.");
    params.addRequiredParam<Real>("permeability_constant", "permeability constant");
    params.addRequiredParam<Real>("solidus_temperature", "Solidus temperature.");
    params.addRequiredParam<Real>("liquidus_temperature", "Liquidus temperature.");
    params.addRequiredParam<Real>("latent_heat", "Latent heat.");
    params.addRequiredParam<Real>("capillary_coefficient", "capillary coefficient.");
    params.addRequiredParam<Real>("thermal_expansion", "thermal expansion.");
    params.addParam<Real>("reference_temperature", 300, "reference_temperature.");
    params.addRequiredParam<Real>("thermalcapillary_coefficient", "thermalcapillary coefficient.");
    params.addParam<RealVectorValue>("gravity", "Direction of the gravity vector");
    // params.addParam<FunctionName>("laser_center_x",
    //                               0,
    //                               "The laser center function of x coordinate.");
    // params.addParam<FunctionName>("laser_center_y",
    //                               0,
    //                               "The laser center function of y coordinate.");
    // params.addParam<FunctionName>("laser_center_z",
    //                               0,
    //                               "The laser center function of z coordinate.");
    params.addParam<UserObjectName>("location", "Location UO");
    params.addRequiredParam<Real>("laser_power", "Laser power.");
    params.addRequiredParam<Real>("effective_beam_radius", "Effective beam radius.");
    params.addRequiredParam<Real>("absorption_coefficient", "Absorption coefficient.");
    params.addRequiredParam<Real>("heat_transfer_coefficient", "Heat transfer coefficient.");
    params.addRequiredParam<Real>("StefanBoltzmann_constant", "Stefan Boltzmann constant.");
    params.addRequiredParam<Real>("material_emissivity", "Material emissivity.");
    params.addRequiredParam<Real>("ambient_temperature", "Ambient temperature.");
    params.addParam<Real>("mass_rate", 2.5e-4, "Mass rate.");
    params.addParam<Real>("mass_radius", 0.25e-3, "Mass radius."););

template <ComputeStage compute_stage>
DEDMaterial<compute_stage>::DEDMaterial(const InputParameters & parameters)
  : ADMaterial<compute_stage>(parameters),
    _velocity(adCoupledVectorValue("velocity")),
    _ls(adCoupledValue("level_set")),
    _temp(adCoupledValue("temperature")),
    _grad_temp(adCoupledGradient("temperature")),
    _grad_ls1(adCoupledGradient("level_set")),
    _grad_ls(adCoupledVectorValue("grad_level_set")),
    _curvature(adCoupledValue("curvature")),
    _rho(declareADProperty<Real>("rho")),
    _grad_rho(declareADProperty<RealVectorValue>("grad_rho")),
    _mu(declareADProperty<Real>("mu")),
    _h(declareADProperty<Real>("enthalpy")),
    _k(declareADProperty<Real>("thermal_conductivity")),
    _dhdT(declareADProperty<Real>("dhdT")),
    _rho_m(declareADProperty<Real>("rho_ixture")),
    _h_m(declareADProperty<Real>("enthalpy_mixture")),
    _mu_m(declareADProperty<Real>("mu_mixture")),
    _powder_feed(declareADProperty<Real>("powder_feed")),
    _mass_change(declareADProperty<Real>("mass_change")),
    _ded_momentum(declareADProperty<RealVectorValue>("ded_momentum")),
    _grads_T(declareADProperty<RealVectorValue>("grads_T")),
    _permeability(declareADProperty<Real>("permeability")),
    _gravity(getParam<RealVectorValue>("gravity")),
    _thermal_expansion(getParam<Real>("thermal_expansion")),
    _reference_temperature(getParam<Real>("reference_temperature")),
    _f_l(declareADProperty<Real>("liquid_mass_fraction")),
    _f_s(declareADProperty<Real>("solid_mass_fraction")),
    _rho_g(getParam<Real>("rho_g")),
    _rho_s(getParam<Real>("rho_s")),
    _rho_l(getParam<Real>("rho_l")),
    _mu_g(getParam<Real>("mu_g")),
    _mu_l(getParam<Real>("mu_l")),
    _mu_s(getParam<Real>("mu_s")),
    _c_g(getParam<Real>("c_g")),
    _c_s(getParam<Real>("c_s")),
    _c_l(getParam<Real>("c_l")),
    _k_g(getParam<Real>("k_g")),
    _k_s(getParam<Real>("k_s")),
    _k_l(getParam<Real>("k_l")),
    _solidus_temperature(getParam<Real>("solidus_temperature")),
    _liquidus_temperature(getParam<Real>("liquidus_temperature")),
    _latent_heat(getParam<Real>("latent_heat")),
    _K0(getParam<Real>("permeability_constant")),
    _sigma(getParam<Real>("capillary_coefficient")),
    _sigmaT(getParam<Real>("thermalcapillary_coefficient")),
    // _laser_center_x(getFunction("laser_center_x")),
    // _laser_center_y(getFunction("laser_center_y")),
    // _laser_center_z(getFunction("laser_center_z")),
    _location(getUserObjectByName<DEDLevelSetLocation>("location")),
    _power(getParam<Real>("laser_power")),
    _alpha(getParam<Real>("absorption_coefficient")),
    _Rb(getParam<Real>("effective_beam_radius")),
    _Ah(getParam<Real>("heat_transfer_coefficient")),
    _stefan_boltzmann(getParam<Real>("StefanBoltzmann_constant")),
    _varepsilon(getParam<Real>("material_emissivity")),
    _T0(getParam<Real>("ambient_temperature")),
    _heat_source(declareADProperty<Real>("heat_source")),
    _mass_rate(getParam<Real>("mass_rate")),
    _mass_radius(getParam<Real>("mass_radius"))
{
}

template <ComputeStage compute_stage>
void
DEDMaterial<compute_stage>::computeQpProperties()
{
  _f_l[_qp] = 1;
  if (_temp[_qp] < _solidus_temperature)
    _f_l[_qp] = 0;
  else if (_temp[_qp] >= _solidus_temperature && _temp[_qp] <= _liquidus_temperature)
    _f_l[_qp] =
        (_temp[_qp] - _solidus_temperature) / (_liquidus_temperature - _solidus_temperature);

  _f_s[_qp] = 1.0 - _f_l[_qp];

  ADReal g_l = _f_l[_qp] / _rho_l / ((1 - _f_l[_qp]) / _rho_s + _f_l[_qp] / _rho_l);
  ADReal g_s = 1 - g_l;

  _f_s[_qp] *= (1 - _ls[_qp]);
  _f_l[_qp] *= (1 - _ls[_qp]);

  // g_l *= (1 - _ls[_qp]);
  // g_s *= (1 - _ls[_qp]);

  ADReal delta_l = (_c_s - _c_l) * _solidus_temperature + _latent_heat;

  ADReal k_s = 0.0138 * _temp[_qp] + 9.13;

  _rho_m[_qp] = g_s * _rho_s + g_l * _rho_l;
  ADReal c_m = _f_s[_qp] * _c_s + _f_l[_qp] * _c_l;
  ADReal k_m = 1.0 / (g_s / k_s + g_l / _k_l);
  _h_m[_qp] = c_m * _temp[_qp] + _f_l[_qp] * delta_l;
  ADReal h_g = _c_g * _temp[_qp];

  //_mu_m[_qp] = _f_s[_qp] * _mu_s + _f_l[_qp] * _mu_l;
  _mu_m[_qp] = _mu_l * _rho_m[_qp] / _rho_l;

  _rho[_qp] = (1 - _ls[_qp]) * _rho_m[_qp] + _ls[_qp] * _rho_g;
  _grad_rho[_qp] = -_grad_ls[_qp] * _rho_m[_qp] + _grad_ls[_qp] * _rho_g;

  _mu[_qp] = (1 - _ls[_qp]) * _mu_m[_qp] + _ls[_qp] * _mu_g;
  _h[_qp] = (1 - _ls[_qp]) * _h_m[_qp] + _ls[_qp] * h_g;
  _k[_qp] = (1 - _ls[_qp]) * k_m + _ls[_qp] * _k_g;
  // if (_temp[_qp] >= _solidus_temperature && _temp[_qp] <= _liquidus_temperature)
  //   _dhdT[_qp] =
  //       (1 - _ls[_qp]) * (c_m + 1.0 / (_liquidus_temperature - _solidus_temperature) * delta_l) +
  //       _ls[_qp] * _c_g;
  // else
  _dhdT[_qp] = (1 - _ls[_qp]) * c_m + _ls[_qp] * _c_g;

  // ADRealVectorValue laser_center(_laser_center_x.value(_t, _q_point[_qp]),
  //                                _laser_center_y.value(_t, _q_point[_qp]),
  //                                _laser_center_z.value(_t, _q_point[_qp]));
  ADRealVectorValue laser_center = _location.getLaserSpotLocation();
  ADReal r = (_q_point[_qp] - laser_center).norm();

  // ADReal heat_term1 =
  //     (_rho[_qp] * (_c_g * _temp[_qp] - _h_m[_qp]) + _h[_qp] * (_rho_g - _rho_m[_qp])) *
  //     _ls_dot[_qp];
  ADReal heat_term2 = 2 * _power * _alpha / (libMesh::pi * Utility::pow<2>(_Rb)) *
                      std::exp(-2.0 * Utility::pow<2>(r / _Rb));
  // if (_t < 0.05)
  //   heat_term2 = 0.0;

  ADReal heat_term3 = -_Ah * (_temp[_qp] - _T0);
  ADReal heat_term4 =
      -_stefan_boltzmann * _varepsilon * (Utility::pow<4>(_temp[_qp]) - Utility::pow<4>(_T0));

  // return -_test[_i][_qp] * (term1 + term2 + term3 + term4) * (_grad_ls[_qp] + epsilon).norm();
  // return -_test[_i][_qp] * (term2) * (_grad_ls[_qp] + epsilon).norm();

  _heat_source[_qp] =
      (heat_term2 + heat_term3 + heat_term4) * (_grad_ls[_qp] + RealVectorValue(1.0e-10)).norm();

  // Real R = 0.1e-3;
  ADReal R = _mass_radius;
  ADReal factor = _mass_rate; // 2.5e-3;
  if (r > R)
    _powder_feed[_qp] = 0.0;
  else
    _powder_feed[_qp] = factor * std::exp(-Utility::pow<2>(r / R));

  _permeability[_qp] = 1 / _K0 * Utility::pow<2>(1 - _f_l[_qp]) /
                       (Utility::pow<3>(_f_l[_qp]) + 1.0e-3) * std::abs(1 - _ls[_qp]);
  // _mass_change[_qp] = _grad_ls[_qp].norm() * (_powder_feed[_qp]);

  _mass_change[_qp] = (_grad_ls[_qp] + RealVectorValue(1.0e-10)).norm() * (_powder_feed[_qp]);

  ADRealVectorValue normal = _grad_ls[_qp] / (_grad_ls[_qp] + RealVectorValue(1e-10)).norm();
  RankTwoTensor iden(RankTwoTensor::initIdentity);
  ADRankTwoTensor proj;
  proj.vectorOuterProduct(normal, normal);
  proj = iden - proj;

  _grads_T[_qp] = proj * _grad_temp[_qp];

  ADRealVectorValue term1 = -_mu_m[_qp] * _permeability[_qp] * _velocity[_qp];
  ADRealVectorValue term2 = _sigma * _curvature[_qp] * _grad_ls[_qp];
  ADRealVectorValue term3 =
      -proj * _grad_temp[_qp] * _sigmaT * (_grad_ls[_qp] + RealVectorValue(1e-10)).norm();
  ADRealVectorValue term4 = -_mass_change[_qp] * _velocity[_qp];
  ADRealVectorValue term5 = -_rho_l * _gravity * _thermal_expansion *
                            (_temp[_qp] - _reference_temperature) * (1 - _ls[_qp]);
  // _ded_momentum[_qp] = term1 + term2 + term3 + term5;
  _ded_momentum[_qp] = -term3 + term2 + term5 + term1;
}
