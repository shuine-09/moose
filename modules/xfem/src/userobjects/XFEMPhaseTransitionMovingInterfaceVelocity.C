//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "XFEMPhaseTransitionMovingInterfaceVelocity.h"

registerMooseObject("XFEMApp", XFEMPhaseTransitionMovingInterfaceVelocity);

InputParameters
XFEMPhaseTransitionMovingInterfaceVelocity::validParams()
{
  InputParameters params = XFEMMovingInterfaceVelocityBase::validParams();
  params.addRequiredParam<Real>("diffusivity_at_positive_level_set",
                                "Diffusivity for level set positive region.");
  params.addRequiredParam<Real>("diffusivity_at_negative_level_set",
                                "Diffusivity for level set negative region.");
  params.addRequiredParam<Real>("equilibrium_concentration_jump",
                                "The jump of the equilibrium concentration at the interface.");
  params.addClassDescription(
      "calculate the interface velocity for a simple phase transition problem.");
  return params;
}

XFEMPhaseTransitionMovingInterfaceVelocity::XFEMPhaseTransitionMovingInterfaceVelocity(
    const InputParameters & parameters)
  : XFEMMovingInterfaceVelocityBase(parameters),
    _diffusivity_at_positive_level_set(getParam<Real>("diffusivity_at_positive_level_set")),
    _diffusivity_at_negative_level_set(getParam<Real>("diffusivity_at_negative_level_set")),
    _equilibrium_concentration_jump(getParam<Real>("equilibrium_concentration_jump"))
{
}

Real
XFEMPhaseTransitionMovingInterfaceVelocity::computeMovingInterfaceVelocity(
    unsigned int point_id) const
{
  Real value_positive = _value_at_interface_uo->getValueAtPositiveLevelSet()[point_id];
  Real value_negative = _value_at_interface_uo->getValueAtNegativeLevelSet()[point_id];
  RealVectorValue grad_positive = _value_at_interface_uo->getGradientAtPositiveLevelSet()[point_id];
  RealVectorValue grad_negative = _value_at_interface_uo->getGradientAtNegativeLevelSet()[point_id];

  // Real xt = (_value_at_interface_uo->getPointCurrentLocation(point_id))(0);
  //
  // std::cout << "xt: " << xt << std::endl;
  std::cout << "value_positive : " << value_positive << std::endl;
  std::cout << "value_negative : " << value_negative << std::endl;
  std::cout << "grad_positive : " << grad_positive(0) << std::endl;
  std::cout << "grad_negative : " << grad_negative(0) << std::endl;

  // Current implementation only supports the case that the interface is moving horizontally
  return std::abs((_diffusivity_at_positive_level_set * grad_positive(0) -
                   _diffusivity_at_negative_level_set * grad_negative(0)) /
                  (value_positive - value_negative + _equilibrium_concentration_jump));
}
