/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
//  Crystal plasticity state variable evolution rate userobject class.
//
#include "CrystalPlasticityStateVariableEvolutionRate.h"

template<>
InputParameters validParams<CrystalPlasticityStateVariableEvolutionRate>()
{
  InputParameters params = validParams<CrystalPlasticityUOBase>();
  params.addRequiredParam<int >("nss", "Number of slip systems");
  params.addParam<std::vector<std::string> >("uo_state_var_evol_rate_comp_name", "Name of state variable evolution rate component property: Same as state variable evolution rate component user object specified in input file.");
  params.addClassDescription("Crystal plasticity state variable evolution rate class.  Override the virtual functions in your class");
  return params;
}

CrystalPlasticityStateVariableEvolutionRate::CrystalPlasticityStateVariableEvolutionRate(const InputParameters & parameters) :
  CrystalPlasticityUOBase(parameters),
  _nss(getParam<int>("nss")),
  _num_mat_state_var_evol_rate_comps(parameters.get<std::vector<std::string> >("uo_state_var_evol_rate_comp_name").size())
{
  _mat_prop_state_var_evol_rate_comps.resize(_num_mat_state_var_evol_rate_comps);

  for (unsigned int i = 0 ; i < _num_mat_state_var_evol_rate_comps; ++i)
  {   
    _mat_prop_state_var_evol_rate_comps[i] = &getMaterialProperty< std::vector<Real> >(parameters.get<std::vector<std::string> >("uo_state_var_evol_rate_comp_name")[i]);
  }
}

bool
CrystalPlasticityStateVariableEvolutionRate::calcStateVariableEvolutionRate(unsigned int qp, std::vector<Real> & val) const
{
  val.resize(_nss);
  for (unsigned int i = 0; i < val.size(); i++)
    val[i] = 0.0;

  for (unsigned int i = 0; i < _num_mat_state_var_evol_rate_comps; i++)
    for (unsigned int j = 0; j < _nss; j++)
      val[j] += (*_mat_prop_state_var_evol_rate_comps[i])[qp][j];

  return true;
}

std::string
CrystalPlasticityStateVariableEvolutionRate::crystalPlasticityUOName() const
{
  return "uo_state_evol_rate_comp_" + _uo_name;
}
