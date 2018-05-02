//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeWeldingIsotropicElasticityTensor.h"
#include "WeldStateIndicator.h"

registerMooseObject("TensorMechanicsApp", ComputeWeldingIsotropicElasticityTensor);

template <>
InputParameters
validParams<ComputeWeldingIsotropicElasticityTensor>()
{
  InputParameters params = validParams<ComputeIsotropicElasticityTensor>();
  params.addClassDescription("Compute weld isotropic elasticity tensor.");
  params.addRequiredParam<UserObjectName>("weld_state", "UserObject name of weld state indicator.");
  return params;
}

ComputeWeldingIsotropicElasticityTensor::ComputeWeldingIsotropicElasticityTensor(
    const InputParameters & parameters)
  : ComputeIsotropicElasticityTensor(parameters)
{
  UserObjectName uo_name = getParam<UserObjectName>("weld_state");
  const UserObject * uo = &(getUserObjectByName<WeldStateIndicator>(uo_name));

  if (dynamic_cast<const WeldStateIndicator *>(uo) == nullptr)
    mooseError(
        "UserObject casting to WeldStateIndicator in ComputeWeldingIsotropicElasticityTensor");

  _weld_state = dynamic_cast<const WeldStateIndicator *>(uo);
}

void
ComputeWeldingIsotropicElasticityTensor::computeQpElasticityTensor()
{
  WeldStateIndicator::WeldStateType _state = _weld_state->getWeldState();

  Real factor = 1.0;
  // Assign elasticity tensor at a given quad point
  if (_state == WeldStateIndicator::WeldStateType::BEFORE)
    factor = 0.01;
  else if (_state == WeldStateIndicator::WeldStateType::HEATING)
    factor = 0.01;
  else if (_state == WeldStateIndicator::WeldStateType::COOLING)
    factor = 1.0;
  else
    factor = 1.0;

  _elasticity_tensor[_qp] = factor * _Cijkl;
}
