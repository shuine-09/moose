/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "CPDislocationDensityRateCompGeneral.h"

template<>
InputParameters validParams<CPDislocationDensityRateCompGeneral>()
{
  InputParameters params = validParams<CrystalPlasticityStateVarRateComponent>();
  params.addRequiredParam<std::string>("uo_dislocation_density_name", "Name of dislocation density property: Same as state variable user object specified in input file.");
  params.addRequiredParam<std::string>("uo_rate_name", "Name of shear rate property: Same as state variable user object specified in input file.");
  params.addRequiredParam<Real>("prefactor", "Prefactor value");
  params.addClassDescription("Dislocation based constitutive model userobject class for dislocation density rate component for a general shear mechanism. Override the virtual functions in your class");
  return params;
}

CPDislocationDensityRateCompGeneral::CPDislocationDensityRateCompGeneral(const InputParameters & parameters) :
    CrystalPlasticityStateVarRateComponent(parameters),
    _mat_prop_dislocation_density(getMaterialProperty<std::vector<Real> >(parameters.get<std::string>("uo_dislocation_density_name"))),
    _mat_prop_rate(getMaterialProperty<std::vector<Real> >(parameters.get<std::string>("uo_rate_name"))),
    _prefactor(getParam<Real>("prefactor"))
{
}

bool
CPDislocationDensityRateCompGeneral::calcStateVariableEvolutionRateComponent(unsigned int qp, std::vector<Real> & val) const
{
  val.assign(_variable_size, 0.0);

  Real tot_rho_m = 0.0;
  for (unsigned int i = 0; i < _variable_size; ++i)
    val[i] = _prefactor * _mat_prop_dislocation_density[qp][i] * std::abs(_mat_prop_rate[qp][i]);

  return true;
}
