//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "XFEMC4VelocityZrOxA.h"
#include "MooseUtils.h"

registerMooseObject("XFEMApp", XFEMC4VelocityZrOxA);

template <>
InputParameters
validParams<XFEMC4VelocityZrOxA>()
{
  InputParameters params = validParams<XFEMMovingInterfaceVelocityBase>();
  params.addParam<Real>("temperature", 1473.15, "Temperature of the cladding (K)");
  params.addClassDescription("Calculate the metal/oxide interface velocity for the 1 interface C4 "
                             "model for Zircaloy-4 corrosion.");
  return params;
}

XFEMC4VelocityZrOxA::XFEMC4VelocityZrOxA(const InputParameters & parameters)
  : XFEMMovingInterfaceVelocityBase(parameters), _temperature(getParam<Real>("temperature"))
{
}

Real
XFEMC4VelocityZrOxA::computeMovingInterfaceVelocity(unsigned int point_id) const
{
  RealVectorValue grad_positive = _value_at_interface_uo->getGradientAtPositiveLevelSet()[point_id];
  RealVectorValue grad_negative = _value_at_interface_uo->getGradientAtNegativeLevelSet()[point_id];

  Real xt = (_value_at_interface_uo->getPointCurrentLocation(0))(0);

  //  std::cout << "point_id: " << point_id << std::endl;

  std::cout << "xt: " << xt << std::endl;

  const Real zirconium_PBR(1.55);
  Real delta = zirconium_PBR * std::abs(xt - 600);
  // std::cout << "delta : " << delta << std::endl;

  // Current implementation only supports the case that the interface is moving horizontally

  // Temperature independent parameters :
  const Real Na(6.022140857e23);
  const Real zircaloy_density(6.56);
  const Real zirconium_atmass(91.22);
  const Real con_zr = zircaloy_density * Na / zirconium_atmass * 1e6;
  const Real zro2_density(5.62);
  const Real zro2_molmass(123.22);
  const Real con_zro2 = zro2_density * Na / zro2_molmass * 1e6;

  const Real Kb(8.6173303e-5);
  const Real migr_jp_f(1e13);
  const Real migr_jp_l(5e-10);
  const Real migr_e_v(1.35);
  const Real migr_e_e(1.30);

  // Temperature dependent parameters
  const Real x_o_ox_m =
      0.667118123 - 1.10606e-5 * _temperature; // the original one, not the weak equivalent

  Real x_o_m_ox = 0.2978;
  if (_temperature > 473.15 && _temperature < 1478.15)
  {
    x_o_m_ox = (28.6 + exp(-6748 / _temperature + 4.748)) * 1e-2;
  }
  else if (_temperature > 1478.15 && _temperature < 1798.15)
  {
    x_o_m_ox = (28.6 + exp(-6301 / _temperature + 4.460)) * 1e-2;
  }
  else if (_temperature > 1798.15 && _temperature < 2338.15)
  {
    x_o_m_ox = (28.6 + exp(-7012 / _temperature + 8.434 - 3.521e-3 * _temperature)) * 1e-2;
  }
  else
  {
    x_o_m_ox = 28.6 * 1e-2;
  }

  const Real con_e_ox_w(1e23);
  const Real con_v_ox_w(1e23);
  const Real con_v_ox_m = con_zro2 * (2 - 3 * x_o_ox_m);
  const Real con_e_ox_m = 2 * con_v_ox_m;

  // Diffusion coefficient is also temperature dependent [um^2/s]
  Real diffusivity_alpha = 10.3;
  if (MooseUtils::absoluteFuzzyEqual(_temperature, 633.15, 1))
  {
    diffusivity_alpha = 1.36e-7;
  }
  else if (MooseUtils::absoluteFuzzyEqual(_temperature, 1223.15, 1))
  {
    diffusivity_alpha = 0.43;
  }
  else if (MooseUtils::absoluteFuzzyEqual(_temperature, 1273.15, 1))
  {
    diffusivity_alpha = 0.3807;
  }
  else if (MooseUtils::absoluteFuzzyEqual(_temperature, 1373.15, 1))
  {
    diffusivity_alpha = 2.40;
  }
  else if (MooseUtils::absoluteFuzzyEqual(_temperature, 1473.15, 1))
  {
    diffusivity_alpha = 10.3;
  }
  else if (MooseUtils::absoluteFuzzyEqual(_temperature, 1573.15, 1))
  {
    diffusivity_alpha = 30;
  }
  else if (MooseUtils::absoluteFuzzyEqual(_temperature, 1673.15, 1))
  {
    diffusivity_alpha = 75.25;
  }
  else if (MooseUtils::absoluteFuzzyEqual(_temperature, 1773.15, 1))
  {
    diffusivity_alpha = 170.25;
  }
  else
  {
    diffusivity_alpha = 9.76 * exp(-54325 / 1.987 / _temperature) * 1e8;
  }

  // Mobilities
  const Real mobil_v = 4 * pow(migr_jp_l, 2) * migr_jp_f * 2 / (Kb * _temperature) *
                       exp(-migr_e_v / (Kb * _temperature));
  const Real mobil_e = 4 * pow(migr_jp_l, 2) * migr_jp_f * -1 / (Kb * _temperature) *
                       exp(-migr_e_e / (Kb * _temperature));

  // Solving for eta
  const Real A = mobil_e * con_e_ox_w - 2 * mobil_v * con_v_ox_m;
  const Real B = mobil_e * con_e_ox_w - mobil_e * con_e_ox_m;
  const Real C = 2 * mobil_v * con_v_ox_w - mobil_e * con_e_ox_m;

  const Real eta = (-B - sqrt(pow(B, 2) - 4 * A * C)) / (2 * A);

  const Real potential = Kb * _temperature * log(eta) / delta;

  _J_v = 1e6 * mobil_v * potential * (con_v_ox_w - con_v_ox_m * pow(eta, 2)) / (1 - pow(eta, 2));

  const Real J_o = -diffusivity_alpha * con_zr * (-grad_negative(0) * 1e-6);

  // std::cout << "ox_a_grad_negative : " << grad_negative(0) << std::endl;
  // std::cout << "ox_a_grad_positive : " << grad_positive(0) << std::endl;

  // std::cout << "J_v : " << J_v << std::endl;
  // std::cout << "J_o : " << J_o << std::endl;

  const Real v_ox_a_init =
      sqrt(0.01126 * exp(-35890 / (1.987 * _temperature)) / (2 * _t)) * (-1e-2);
  const Real v_ox_a = 1e6 * (J_o - _J_v) / (con_zr * (3 * x_o_ox_m - (x_o_m_ox / (1 - x_o_m_ox))));
  //  _ox_vel = v_ox_a;

  if (delta == 0)
  {
    std::cout << "Oxide-alpha (initial) velocity : " << v_ox_a_init << std::endl;
    return v_ox_a_init;
  }
  else
    std::cout << "Oxide-alpha velocity : " << v_ox_a << std::endl;
  return v_ox_a;
}
