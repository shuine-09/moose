//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "XFEMC4VelocityZrABForOpti.h"
#include "MooseUtils.h"

registerMooseObject("XFEMApp", XFEMC4VelocityZrABForOpti);

template <>
InputParameters
validParams<XFEMC4VelocityZrABForOpti>()
{
  InputParameters params = validParams<XFEMMovingInterfaceVelocityBase>();
  params.addParam<Real>("diffusivity_alpha",10,"oxygen diffusivity in the alpha phase (µm²/s)");
  params.addParam<Real>("diffusivity_beta",60,"oxygen diffusivity in the beta phase (µm²/s)");
  params.addParam<Real>("temperature", 1473.15, "Temperature of the cladding (K)");
  params.addClassDescription(
      "Calculate the alpha phase/beta phase interface velocity for the 2 interfaces C4 model for Zircaloy-4 corrosion.");
  return params;
}

XFEMC4VelocityZrABForOpti::XFEMC4VelocityZrABForOpti(const InputParameters & parameters)
  : XFEMMovingInterfaceVelocityBase(parameters),
    _diffusivity_alpha(getParam<Real>("diffusivity_alpha")),
    _diffusivity_beta(getParam<Real>("diffusivity_beta")),
    _temperature(getParam<Real>("temperature"))
{
}

Real
XFEMC4VelocityZrABForOpti::computeMovingInterfaceVelocity(unsigned int point_id) const
{
  RealVectorValue grad_positive = _value_at_interface_uo->getGradientAtPositiveLevelSet()[point_id];
  RealVectorValue grad_negative = _value_at_interface_uo->getGradientAtNegativeLevelSet()[point_id];

//  Real xt = (_value_at_interface_uo->getPointCurrentLocation(point_id))(0);

  //  std::cout << "point_id: " << point_id << std::endl;

//  std::cout << "xt: " << xt << std::endl;

  // Current implementation only supports the case that the interface is moving horizontally

//  Values at the interface (strong discontinuity)
  Real x_o_b_a = (9.59e-3 * (_temperature - 1136) + 4.72e-6 * pow(_temperature - 1136,2) - 4.35e-9 * pow(_temperature - 1136,3)) * 1e-2;
  Real x_o_a_b = (45.86e-3 * (_temperature - 1136) - 44.77e-6 * pow(_temperature - 1136,2) + 17.40e-9 * pow(_temperature - 1136,3)) * 1e-2; //the original one, not the weak equivalent
  const Real c_o_a_b = x_o_a_b / (1 - x_o_a_b);
  const Real c_o_b_a = x_o_b_a / (1 - x_o_b_a);

  const Real J_b_to_a = -_diffusivity_alpha * grad_positive(0);
  const Real J_a_to_b = -_diffusivity_beta * (-grad_negative(0));

  //std::cout << "Da : " << _diffusivity_alpha << std::endl;
  //std::cout << "Db : " << _diffusivity_beta << std::endl;

  //std::cout << "ab_grad_negative : " << grad_negative(0) << std::endl;
  //std::cout << "ab_grad_positive : " << grad_positive(0) << std::endl;

  //std::cout << "J_b_to_a : " << J_b_to_a * 4.33e28<< std::endl;
  //std::cout << "J_a_to_b : " << J_a_to_b * 4.33e28<< std::endl;

  const Real v_a_b = (J_b_to_a + J_a_to_b) / (c_o_a_b - c_o_b_a);

  //std::cout << "Alpha-beta velocity : " << v_a_b << std::endl;

  return v_a_b;


}
