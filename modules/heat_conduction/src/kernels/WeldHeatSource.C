//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "WeldHeatSource.h"
#include "Function.h"
#include "WeldStateIndicator.h"

registerMooseObject("HeatConductionApp", WeldHeatSource);

template <>
InputParameters
validParams<WeldHeatSource>()
{
  InputParameters params = validParams<BodyForce>();

  // Override defaults and documentation, weak form is identical to BodyForce in MOOSE
  // params.addParam<Real>("value", 1.0, "Value of heat source. Multiplied by function if
  // present.");
  params.addParam<FunctionName>("function", "1", "Function describing the volumetric heat source");
  params.addRequiredParam<UserObjectName>("weld_state", "UserObject name of weld state indicator.");
  return params;
}

WeldHeatSource::WeldHeatSource(const InputParameters & parameters) : BodyForce(parameters)
{
  UserObjectName uo_name = getParam<UserObjectName>("weld_state");
  const UserObject * uo = &(getUserObjectByName<WeldStateIndicator>(uo_name));

  if (dynamic_cast<const WeldStateIndicator *>(uo) == nullptr)
    mooseError(
        "UserObject casting to WeldStateIndicator in ComputeWeldingIsotropicElasticityTensor");

  _weld_state = dynamic_cast<const WeldStateIndicator *>(uo);
}

Real
WeldHeatSource::computeQpResidual()
{
  Real factor = _scale * _postprocessor * _function.value(_t, _q_point[_qp]);

  WeldStateIndicator::WeldStateType _state = _weld_state->getWeldState();

  // Assign elasticity tensor at a given quad point
  if (_state != WeldStateIndicator::WeldStateType::HEATING)
    factor = 0.0;

  // std::cout << "WeldHeatSource : State = " << _state << std::endl;

  return _test[_i][_qp] * -factor;
}
