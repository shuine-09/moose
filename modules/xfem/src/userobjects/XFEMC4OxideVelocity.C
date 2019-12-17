//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "XFEMC4OxideVelocity.h"

registerMooseObject("XFEMApp", XFEMC4OxideVelocity);

template <>
InputParameters
validParams<XFEMC4OxideVelocity>()
{
  InputParameters params = validParams<XFEMMovingInterfaceVelocityBase>();
  params.addRequiredParam<Real>("diffusivity_at_positive_level_set",
                                "Diffusivity for level set positive region.");
  params.addRequiredParam<Real>("diffusivity_at_negative_level_set",
                                "Diffusivity for level set negative region.");
  params.addRequiredParam<Real>("equilibrium_concentration_jump",
                                "The jump of the equilibrium concentration at the interface.");
  params.addParam<Real>("x0", 0, "initial x location.");
  params.addClassDescription(
      "calculate the interface velocity for a simple phase transition problem.");
  return params;
}

XFEMC4OxideVelocity::XFEMC4OxideVelocity(const InputParameters & parameters)
  : XFEMMovingInterfaceVelocityBase(parameters),
    _diffusivity_at_positive_level_set(getParam<Real>("diffusivity_at_positive_level_set")),
    _diffusivity_at_negative_level_set(getParam<Real>("diffusivity_at_negative_level_set")),
    _equilibrium_concentration_jump(getParam<Real>("equilibrium_concentration_jump")),
    _x0(getParam<Real>("x0"))
{
}

Real
XFEMC4OxideVelocity::computeMovingInterfaceVelocity(unsigned int point_id) const
{
  Real value_positive = _value_at_interface_uo->getValueAtPositiveLevelSet()[point_id];
  Real value_negative = _value_at_interface_uo->getValueAtNegativeLevelSet()[point_id];
  RealVectorValue grad_positive = _value_at_interface_uo->getGradientAtPositiveLevelSet()[point_id];
  RealVectorValue grad_negative = _value_at_interface_uo->getGradientAtNegativeLevelSet()[point_id];

  Real xt = (_value_at_interface_uo->getPointCurrentLocation(point_id))(0);

  Real delta = std::abs(xt - _x0);

  //  std::cout << "delta: " << delta << std::endl;

  // Current implementation only supports the case that the interface is moving horizontally
  //  return std::abs((_diffusivity_at_positive_level_set * grad_positive(0) -
  //                   _diffusivity_at_negative_level_set * grad_negative(0)) /
  //                  (value_positive - value_negative + _equilibrium_concentration_jump));
  const Real zircaloy_density(6550);
  const Real zro2_density(5680);
  const Real oxygen_atmass(16);
  const Real zirconium_atmass(91.2);
  const Real zirconium_PBR(1.56);
  const Real Na(6.022140857e23);
  const Real Kb(8.6173303e-5);
  const Real migr_jp_f(1e13);
  const Real migr_jp_l(5e-10);
  const Real migr_e_v(1.35);
  const Real migr_e_e(1.30);
  const Real con_v_ox_w(1e23);
  const Real con_e_ox_w(1e23);
  const Real con_v_ox_m(1.3194e27);
  const Real con_e_ox_m = 2 * con_v_ox_m;
  const Real con_o_ox_m(5.24e28);
  const Real con_o_m_ox(1.7078e28);

  // use fixed temperature at 1200C until we add a temperature diffusion kernel
  const Real temperature(1473);

  const Real mobil_v = 4 * pow(migr_jp_l, 2) * migr_jp_f * 2 / (Kb * temperature) *
                       exp(-migr_e_v / (Kb * temperature));
  const Real mobil_e = 4 * pow(migr_jp_l, 2) * migr_jp_f * -1 / (Kb * temperature) *
                       exp(-migr_e_e / (Kb * temperature));

  const Real A = mobil_e * con_e_ox_w - 2 * mobil_v * con_v_ox_m;
  const Real B = mobil_e * con_e_ox_w - mobil_e * con_e_ox_m;
  const Real C = 2 * mobil_v * con_v_ox_w - mobil_e * con_e_ox_m;
  const Real eta = (-B - sqrt(pow(B, 2) - 4 * A * C)) / (2 * A);
  const Real potential = Kb * temperature * log(eta);

  const Real J_v =
      mobil_v * potential * (con_v_ox_w - con_v_ox_m * pow(eta, 2)) / (1 - pow(eta, 2)) / delta;
  // const Real J_o = zircaloy_density * Na / (zirconium_atmass * 1e-3) *
  //                  _diffusivity_at_positive_level_set * (grad_positive(0) + grad_positive(1) +
  //                  grad_positive(2)) / 3 * 1e-6;
  const Real J_o = zircaloy_density * Na / (zirconium_atmass * 1e-3) *
                   _diffusivity_at_positive_level_set * (grad_positive(0)) * 1e-6;

  if (delta == 0)
    return sqrt(0.01126 * exp(-35890 / (1.987 * temperature)) / (2 * _t)) * (-1e-2);
  else
    return -zirconium_PBR * (J_v - J_o) / (zirconium_PBR * con_o_ox_m - con_o_m_ox);
}
