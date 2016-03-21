/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "CPDislocationDensityRatePptResistance.h"

template<>
InputParameters validParams<CPDislocationDensityRatePptResistance>()
{
  InputParameters params = validParams<CrystalPlasticityStateVarRateComponent>();
  params.addParam<std::string>("uo_rate_name", "Name of shear rate property: Same as state variable user object specified in input file.");
  params.addParam<std::string>("property_in_function_name", "Property name used in function to evaluate property");
  params.addParam<Real>("number_density", 0.0, "Average number density of precipitate");
  params.addParam<Real>("size", 0.0, "Average size of precipitate");
  params.addParam<FunctionName>("factor_function", "Function to obtain shear rate dependent factor");
  return params;
}

CPDislocationDensityRatePptResistance::CPDislocationDensityRatePptResistance(const InputParameters & parameters) :
    CrystalPlasticityStateVarRateComponent(parameters),
    _mat_prop_rate(getMaterialProperty<std::vector<Real> >(parameters.get<std::string>("uo_rate_name"))),
    _mat_prop(getMaterialProperty<std::vector<Real> >(parameters.get<std::string>("property_in_function_name"))),
    _number_density(getParam<Real>("number_density")),
    _size(getParam<Real>("size")),
    _factor_function(isParamValid("factor_function") ? &getFunction("factor_function") : NULL)
{
}

bool
CPDislocationDensityRatePptResistance::calcStateVariableEvolutionRateComponent(unsigned int qp, std::vector<Real> & val) const
{
  Point p;
  for (unsigned int i = 0; i < _variable_size; ++i)
    val[i] = _factor_function->value(_mat_prop[qp][i], p) * std::sqrt(_number_density * _size) * std::abs(_mat_prop_rate[qp][i]);

  return true;
}
