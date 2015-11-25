/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
//  Crystal plasticity state variable evolution rate component userobject class.
//
#include "CrystalPlasticityStateVariableEvolutionRateComponentGSS.h"

template<>
InputParameters validParams<CrystalPlasticityStateVariableEvolutionRateComponentGSS>()
{
  InputParameters params = validParams<CrystalPlasticityStateVariableEvolutionRateComponent>();
  params.addParam<std::vector<std::string> >("uo_slip_rate_name", "Name of slip rate property: Same as slip rate user object specified in input file.");
  params.addClassDescription("Crystal plasticity state variable evolution rate component base class.  Override the virtual functions in your class");
  return params;
}

CrystalPlasticityStateVariableEvolutionRateComponentGSS::CrystalPlasticityStateVariableEvolutionRateComponentGSS(const InputParameters & parameters) :
  CrystalPlasticityStateVariableEvolutionRateComponent(parameters),
  _num_mat_slip_rates(parameters.get<std::vector<std::string> >("uo_slip_rate_name").size())
{
  _mat_prop_slip_rates.resize(_num_mat_slip_rates);

  for (unsigned int i = 0 ; i < _num_mat_slip_rates; ++i)
  {
    _mat_prop_slip_rates[i] = &getMaterialProperty< std::vector<Real> >(parameters.get<std::vector<std::string> >("uo_slip_rate_name")[i]);
  }
}

bool
CrystalPlasticityStateVariableEvolutionRateComponentGSS::calcStateVariableEvolutionRateComponent(unsigned int qp, std::vector<Real> & val) const
{
  val.resize(_nss);
  for (unsigned int i = 0; i < _nss; i++)
    val[i] = (*_mat_prop_slip_rates[0])[qp][i];

  return true;
}
