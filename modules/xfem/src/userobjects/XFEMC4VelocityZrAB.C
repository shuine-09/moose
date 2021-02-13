//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "XFEMC4VelocityZrAB.h"
#include "MooseUtils.h"

registerMooseObject("XFEMApp", XFEMC4VelocityZrAB);

template <>
InputParameters
validParams<XFEMC4VelocityZrAB>()
{
  InputParameters params = validParams<XFEMMovingInterfaceVelocityBase>();
  params.addParam<Real>("temperature", 1473.15, "Temperature of the cladding (K)");
  params.addClassDescription("Calculate the alpha phase/beta phase interface velocity for the 2 "
                             "interfaces C4 model for Zircaloy-4 corrosion.");
  return params;
}

XFEMC4VelocityZrAB::XFEMC4VelocityZrAB(const InputParameters & parameters)
  : XFEMMovingInterfaceVelocityBase(parameters), _temperature(getParam<Real>("temperature"))
{
}

Real
XFEMC4VelocityZrAB::computeMovingInterfaceVelocity(unsigned int point_id) const
{
  RealVectorValue grad_positive = _value_at_interface_uo->getGradientAtPositiveLevelSet()[point_id];
  RealVectorValue grad_negative = _value_at_interface_uo->getGradientAtNegativeLevelSet()[point_id];

  //  Real xt = (_value_at_interface_uo->getPointCurrentLocation(point_id))(0);

  //  std::cout << "point_id: " << point_id << std::endl;

  //  std::cout << "xt: " << xt << std::endl;

  // Current implementation only supports the case that the interface is moving horizontally

  //  Values at the interface (strong discontinuity)
  Real x_o_b_a = (9.59e-3 * (_temperature - 1136) + 4.72e-6 * pow(_temperature - 1136, 2) -
                  4.35e-9 * pow(_temperature - 1136, 3)) *
                 1e-2;
  Real x_o_a_b = (45.86e-3 * (_temperature - 1136) - 44.77e-6 * pow(_temperature - 1136, 2) +
                  17.40e-9 * pow(_temperature - 1136, 3)) *
                 1e-2; // the original one, not the weak equivalent
  const Real c_o_a_b = x_o_a_b / (1 - x_o_a_b);
  const Real c_o_b_a = x_o_b_a / (1 - x_o_b_a);

  // Diffusion coefficients
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
    diffusivity_alpha = 30.0;
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

  Real diffusivity_beta = 280.8594;
  if (MooseUtils::absoluteFuzzyEqual(_temperature, 633.15, 1))
  {
    diffusivity_beta = 4.8e-4;
  }
  else if (MooseUtils::absoluteFuzzyEqual(_temperature, 1223.15, 1))
  {
    diffusivity_beta = 24;
  }
  else if (MooseUtils::absoluteFuzzyEqual(_temperature, 1273.15, 1))
  {
    diffusivity_beta = 30.1;
  }
  else if (MooseUtils::absoluteFuzzyEqual(_temperature, 1373.15, 1))
  {
    diffusivity_beta = 88.8;
  }
  else if (MooseUtils::absoluteFuzzyEqual(_temperature, 1473.15, 1))
  {
    diffusivity_beta = 280.8594;
    // diffusivity_beta = 140;
  }
  else if (MooseUtils::absoluteFuzzyEqual(_temperature, 1573.15, 1))
  {
    diffusivity_beta = 704.6875;
  }
  else if (MooseUtils::absoluteFuzzyEqual(_temperature, 1673.15, 1))
  {
    diffusivity_beta = 709.65625;
  }
  else if (MooseUtils::absoluteFuzzyEqual(_temperature, 1773.15, 1))
  {
    diffusivity_beta = 1360;
  }
  else
  {
    diffusivity_beta = 0.265 * exp(-34205 / 1.987 / _temperature) * 1e8;
  }

  const Real J_b_to_a = -diffusivity_alpha * grad_positive(0);
  const Real J_a_to_b = -diffusivity_beta * (-grad_negative(0));

  std::cout << "ab_grad_negative : " << grad_negative(0) << std::endl;
  std::cout << "ab_grad_positive : " << grad_positive(0) << std::endl;

  // std::cout << "J_a : " << J_b_to_a * 4.33e28<< std::endl;
  // std::cout << "J_b : " << J_a_to_b * 4.33e28<< std::endl;

  const Real v_a_b = (J_b_to_a + J_a_to_b) / (c_o_a_b - c_o_b_a);

  std::cout << "Alpha-beta velocity : " << v_a_b << std::endl;

  return v_a_b;
}
