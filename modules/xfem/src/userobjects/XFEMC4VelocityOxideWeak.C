//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "XFEMC4VelocityOxideWeak.h"

registerMooseObject("XFEMApp", XFEMC4VelocityOxideWeak);

template <>
InputParameters
validParams<XFEMC4VelocityOxideWeak>()
{
  InputParameters params = validParams<XFEMMovingInterfaceVelocityBase>();
  params.addRequiredParam<Real>("diffusivity_alpha",
                                "Diffusivity of oxygen in the Zr alpha phase.");
  params.addParam<Real>("x0", 0, "Initial position of the interface.");
  params.addClassDescription(
      "Calculate the metal/oxide interface velocity for the 1 interface C4 model for Zircaloy-4 corrosion.");
  return params;
}

XFEMC4VelocityOxideWeak::XFEMC4VelocityOxideWeak(const InputParameters & parameters)
  : XFEMMovingInterfaceVelocityBase(parameters),
    _diffusivity_alpha(getParam<Real>("diffusivity_alpha")),
    _x0(getParam<Real>("x0"))
{
}

Real
XFEMC4VelocityOxideWeak::computeMovingInterfaceVelocity(unsigned int point_id) const
{
  RealVectorValue grad_positive = _value_at_interface_uo->getGradientAtPositiveLevelSet()[point_id];
  RealVectorValue grad_negative = _value_at_interface_uo->getGradientAtNegativeLevelSet()[point_id];

  Real xt = (_value_at_interface_uo->getPointCurrentLocation(point_id))(0);

//  std::cout << "xt: " << xt << std::endl;

  Real delta = std::abs(xt - 6e-4);
  std::cout << "delta : " << delta << std::endl;


  // Current implementation only supports the case that the interface is moving horizontally

  const Real Na(6.022140857e23);
  const Real zircaloy_density(6.56);
  const Real zirconium_atmass(91.22);
  const Real con_zr = zircaloy_density * Na / zirconium_atmass * 1e6;
  const Real zro2_density(5.62);
  const Real zro2_molmass(123.22);
  const Real con_zro2 = zro2_density * Na / zro2_molmass * 1e6;
  //const Real zirconium_PBR(1.56);

  const Real Kb(8.6173303e-5);
  const Real migr_jp_f(1e13);
  const Real migr_jp_l(5e-10);
  const Real migr_e_v(1.35);
  const Real migr_e_e(1.30);

  const Real x_o_m_ox(0.2978); //retrieve from constraint ?
  const Real x_o_ox_m(0.6508); //the original one, not the weak equivalent
  //const Real x_o_ox(0.6667);  //the original one, not the weak equivalent

  const Real con_e_ox_w(1e23);
  const Real con_v_ox_w(1e23);
  const Real con_v_ox_m = con_zro2 * (2 - 3 * x_o_ox_m);
  const Real con_e_ox_m = 2 * con_v_ox_m;

// use fixed temperature at 1200C until we add a temperature diffusion kernel
  const Real temperature(1473);

  const Real mobil_v = 4 * pow(migr_jp_l, 2) * migr_jp_f * 2 / (Kb * temperature) *
                       exp(-migr_e_v / (Kb * temperature));
  const Real mobil_e = 4 * pow(migr_jp_l, 2) * migr_jp_f * -1 / (Kb * temperature) *
                       exp(-migr_e_e / (Kb * temperature));

  const Real A = mobil_e * con_e_ox_w - 2 * mobil_v * con_v_ox_m;
  const Real B = mobil_e * con_e_ox_w - mobil_e * con_e_ox_m;
  const Real C = 2 * mobil_v * con_v_ox_w - mobil_e * con_e_ox_m;

  const Real eta = (- B - sqrt(pow(B, 2) - 4 * A * C)) / (2 * A);

  const Real potential = Kb * temperature * log(eta)/delta;

  const Real J_v = mobil_v * potential * (con_v_ox_w - con_v_ox_m * pow(eta,2)) / (1 - pow(eta,2));
  const Real J_o = -_diffusivity_alpha * con_zr / pow(1 - x_o_m_ox, 2) * (-grad_negative(0));

  std::cout << "J_v : " << J_v << std::endl;
  std::cout << "J_o : " << J_o << std::endl;

  const Real v_ox_a_init = sqrt(0.01126 * exp(-35890 / (1.987 * temperature)) / (2 * _t)) * (-1e-2);
  const Real v_ox_a = (J_o - J_v) / (con_zr * (3 * x_o_ox_m - (x_o_m_ox / (1 - x_o_m_ox))));

  if (delta == 0)
    {
    std::cout << "Oxide-alpha (initial) velocity : " << v_ox_a_init << std::endl;
    return v_ox_a_init;
    }
  else
    std::cout << "Oxide-alpha velocity : " << v_ox_a << std::endl;
    return v_ox_a;


}
