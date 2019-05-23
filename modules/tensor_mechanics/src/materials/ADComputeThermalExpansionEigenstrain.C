//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADComputeThermalExpansionEigenstrain.h"

registerMooseObject("TensorMechanicsApp", ADComputeThermalExpansionEigenstrain);

InputParameters
ADComputeThermalExpansionEigenstrain::validParams()
{
  InputParameters params = ADComputeThermalExpansionEigenstrainBase::validParams();
  params.addClassDescription("Computes eigenstrain due to thermal expansion "
                             "with a constant coefficient");
  params.addRequiredParam<Real>("thermal_expansion_coeff", "Thermal expansion coefficient");
  params.addCoupledVar(
      "activated_elem_aux", 1, "Temperature aux variable used to determine activated elements.");
  return params;
}

ADComputeThermalExpansionEigenstrain::ADComputeThermalExpansionEigenstrain(
    const InputParameters & parameters)
  : ADComputeThermalExpansionEigenstrainBase(parameters),
    _thermal_expansion_coeff(getParam<Real>("thermal_expansion_coeff")),
    _activated_elem(coupledValue("activated_elem_aux"))
{
}

void
ADComputeThermalExpansionEigenstrain::computeThermalStrain(ADReal & thermal_strain)
{
  if (_activated_elem[_qp] >= 1.0 && _current_elem->subdomain_id() == 1)
    thermal_strain =
        _thermal_expansion_coeff * (_temperature[_qp] - _stress_free_temperature[_qp] - 600);
  else
    thermal_strain = _thermal_expansion_coeff * (_temperature[_qp] - _stress_free_temperature[_qp]);
}
