/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "CPDislocationDensityRateCompGeneral.h"

registerMooseObject("TensorMechanicsApp", CPDislocationDensityRateCompGeneral);

template <>
InputParameters
validParams<CPDislocationDensityRateCompGeneral>()
{
  InputParameters params = validParams<CrystalPlasticityStateVarRateComponent>();
  params.addParam<std::string>("uo_prop_name",
                               "Name of dislocation density property: Same as state variable user "
                               "object specified in input file.");
  params.addParam<std::string>(
      "uo_rate_name",
      "Name of shear rate property: Same as state variable user object specified in input file.");
  params.addParam<Real>("prefactor", "Prefactor value");
  params.addParam<bool>("use_dislocation_density", false, "Use dislocation density term.");
  params.addClassDescription(
      "Dislocation based constitutive model userobject class for dislocation density rate "
      "component for a general shear mechanism. Override the virtual functions in your class");
  return params;
}

CPDislocationDensityRateCompGeneral::CPDislocationDensityRateCompGeneral(
    const InputParameters & parameters)
  : CrystalPlasticityStateVarRateComponent(parameters),
    _mat_prop(getMaterialProperty<std::vector<Real>>(parameters.get<std::string>("uo_prop_name"))),
    _mat_prop_rate(
        getMaterialProperty<std::vector<Real>>(parameters.get<std::string>("uo_rate_name"))),
    _prefactor(getParam<Real>("prefactor")),
    _use_dislocation_density(getParam<bool>("use_dislocation_density"))
{
}

bool
CPDislocationDensityRateCompGeneral::calcStateVariableEvolutionRateComponent(
    unsigned int qp, std::vector<Real> & val) const
{
  val.assign(_variable_size, 0.0);

  for (unsigned int i = 0; i < _variable_size; ++i)
    if (_use_dislocation_density)
      val[i] = _prefactor * _mat_prop[qp][i] * std::abs(_mat_prop_rate[qp][i]);
    else
      val[i] = _prefactor * std::abs(_mat_prop_rate[qp][i]);

  return true;
}
