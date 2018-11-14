//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "WeldingHeatConductionMaterial.h"
#include "Function.h"
#include "WeldStateIndicator.h"

#include "libmesh/quadrature.h"

registerMooseObject("HeatConductionApp", WeldingHeatConductionMaterial);

template <>
InputParameters
validParams<WeldingHeatConductionMaterial>()
{
  InputParameters params = validParams<HeatConductionMaterial>();
  params.addRequiredParam<UserObjectName>("weld_state", "UserObject name of weld state indicator.");
  return params;
}

WeldingHeatConductionMaterial::WeldingHeatConductionMaterial(const InputParameters & parameters)
  : HeatConductionMaterial(parameters)
{
  UserObjectName uo_name = getParam<UserObjectName>("weld_state");
  const UserObject * uo = &(getUserObjectByName<WeldStateIndicator>(uo_name));

  if (dynamic_cast<const WeldStateIndicator *>(uo) == nullptr)
    mooseError(
        "UserObject casting to WeldStateIndicator in ComputeWeldingIsotropicElasticityTensor");

  _weld_state = dynamic_cast<const WeldStateIndicator *>(uo);
}

void
WeldingHeatConductionMaterial::computeProperties()
{
  WeldStateIndicator::WeldStateType _state = _weld_state->getWeldState();

  Real factor = 1.0;
  // Assign elasticity tensor at a given quad point
  if (_state == WeldStateIndicator::WeldStateType::BEFORE)
    factor = 0.00001;
  else if (_state == WeldStateIndicator::WeldStateType::HEATING)
    factor = 1.0;
  else if (_state == WeldStateIndicator::WeldStateType::COOLING)
    factor = 1.0;
  else
    factor = 1.0;

  for (unsigned int qp(0); qp < _qrule->n_points(); ++qp)
  {
    Real qp_temperature = 0;
    if (_has_temp)
    {
      qp_temperature = _temperature[qp];
      if (_temperature[qp] < 0)
      {
        std::stringstream msg;
        msg << "WARNING:  In WeldingHeatConductionMaterial:  negative temperature!\n"
            << "\tResetting to zero.\n"
            << "\t_qp: " << qp << "\n"
            << "\ttemp: " << _temperature[qp] << "\n"
            << "\telem: " << _current_elem->id() << "\n"
            << "\tproc: " << processor_id() << "\n";
        mooseWarning(msg.str());
        qp_temperature = 0;
      }
    }
    if (_thermal_conductivity_temperature_function)
    {
      Point p;
      _thermal_conductivity[qp] =
          factor * _thermal_conductivity_temperature_function->value(qp_temperature, p);
      _thermal_conductivity_dT[qp] =
          factor * _thermal_conductivity_temperature_function->timeDerivative(qp_temperature, p);
    }
    else
    {
      _thermal_conductivity[qp] = factor * _my_thermal_conductivity;
      _thermal_conductivity_dT[qp] = 0;
    }

    if (_specific_heat_temperature_function)
    {
      Point p;
      _specific_heat[qp] = factor * _specific_heat_temperature_function->value(qp_temperature, p);
    }
    else
    {
      _specific_heat[qp] = factor * _my_specific_heat;
    }
  }
}
