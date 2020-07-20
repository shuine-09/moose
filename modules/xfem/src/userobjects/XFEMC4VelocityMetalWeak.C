//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "XFEMC4VelocityMetalWeak.h"

registerMooseObject("XFEMApp", XFEMC4VelocityMetalWeak);

template <>
InputParameters
validParams<XFEMC4VelocityMetalWeak>()
{
  InputParameters params = validParams<XFEMMovingInterfaceVelocityBase>();
  params.addRequiredParam<Real>("diffusivity_alpha",
                                "Diffusivity of oxygen in the alpha phase.");
  params.addRequiredParam<Real>("diffusivity_beta",
                                "Diffusivity of oxygen in the beta phase.");
  //params.addParam<Real>("x0", 0, "Initial x location.");
  params.addClassDescription(
      "Calculate the alpha phase/beta phase interface velocity for the 2 interfaces C4 model for Zircaloy-4 corrosion.");
  return params;
}

XFEMC4VelocityMetalWeak::XFEMC4VelocityMetalWeak(const InputParameters & parameters)
  : XFEMMovingInterfaceVelocityBase(parameters),
    _diffusivity_alpha(getParam<Real>("diffusivity_alpha")),
    _diffusivity_beta(getParam<Real>("diffusivity_beta"))
    //_x0(getParam<Real>("x0"))
{
}

Real
XFEMC4VelocityMetalWeak::computeMovingInterfaceVelocity(unsigned int point_id) const
{
  RealVectorValue grad_positive = _value_at_interface_uo->getGradientAtPositiveLevelSet()[point_id];
  RealVectorValue grad_negative = _value_at_interface_uo->getGradientAtNegativeLevelSet()[point_id];

  Real xt = (_value_at_interface_uo->getPointCurrentLocation(point_id))(0);

//  std::cout << "xt: " << xt << std::endl;

  // Current implementation only supports the case that the interface is moving horizontally

  const Real x_o_b_a(0.0360); //retrieve from constraint ?
  const Real x_o_a_b(0.1104); //the original one, not the weak equivalent

// use fixed temperature at 1200C until we add a temperature diffusion kernel
  //const Real temperature(1473);


  const Real J_b_to_a = -_diffusivity_alpha * grad_positive(0);
  const Real J_a_to_b = -_diffusivity_beta * (-grad_negative(0));

  const Real c_o_a_b = x_o_a_b / (1 - x_o_a_b);
  const Real c_o_b_a = x_o_b_a / (1 - x_o_b_a);

  std::cout << "ab_grad_negative : " << grad_negative(0) << std::endl;
  std::cout << "ab_grad_positive : " << grad_positive(0) << std::endl;

  //std::cout << "J_b_to_a : " << J_b_to_a * 4.33e28<< std::endl;
  //std::cout << "J_a_to_b : " << J_a_to_b * 4.33e28<< std::endl;

  const Real v_a_b = (J_b_to_a - J_a_to_b) / (c_o_a_b - c_o_b_a);

  //std::cout << "Alpha-beta velocity : " << v_a_b << std::endl;

  return v_a_b;


}
