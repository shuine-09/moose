//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeWeldThermalExpansionEigenstrain.h"
#include "Function.h"
#include "WeldStateIndicator.h"

registerMooseObject("TensorMechanicsApp", ComputeWeldThermalExpansionEigenstrain);

template <>
InputParameters
validParams<ComputeWeldThermalExpansionEigenstrain>()
{
  InputParameters params = validParams<ComputeThermalExpansionEigenstrainBase>();
  params.addClassDescription("Computes eigenstrain due to thermal expansion "
                             "with a constant coefficient");
  params.addRequiredParam<Real>("thermal_expansion_coeff", "Thermal expansion coefficient");
  params.addRequiredParam<UserObjectName>("weld_state", "UserObject name of weld state indicator.");

  return params;
}

ComputeWeldThermalExpansionEigenstrain::ComputeWeldThermalExpansionEigenstrain(
    const InputParameters & parameters)
  : ComputeThermalExpansionEigenstrainBase(parameters),
    _thermal_expansion_coeff(getParam<Real>("thermal_expansion_coeff"))
{
  UserObjectName uo_name = getParam<UserObjectName>("weld_state");
  const UserObject * uo = &(getUserObjectByName<WeldStateIndicator>(uo_name));

  if (dynamic_cast<const WeldStateIndicator *>(uo) == nullptr)
    mooseError(
        "UserObject casting to WeldStateIndicator in ComputeWeldingIsotropicElasticityTensor");

  _weld_state = dynamic_cast<const WeldStateIndicator *>(uo);
}

void
ComputeWeldThermalExpansionEigenstrain::computeThermalStrain(Real & thermal_strain,
                                                             Real & instantaneous_cte)
{
  WeldStateIndicator::WeldStateType _state = _weld_state->getWeldState();

  Real factor = 1.0;
  // Assign elasticity tensor at a given quad point
  if (_state == WeldStateIndicator::WeldStateType::BEFORE)
  {
    factor = 0.00000001;
  }
  else if (_state == WeldStateIndicator::WeldStateType::HEATING)
  {
    factor = 0.00000001;
  }
  else if (_state == WeldStateIndicator::WeldStateType::COOLING)
  {
    factor = 1.0;
  }
  else
    factor = 1.0;

  thermal_strain =
      _thermal_expansion_coeff * (_temperature[_qp] - _stress_free_temperature[_qp]) * factor;
  instantaneous_cte = _thermal_expansion_coeff * factor;
}
