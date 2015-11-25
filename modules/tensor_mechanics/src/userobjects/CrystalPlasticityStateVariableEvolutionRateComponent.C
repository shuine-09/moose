/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
//  Crystal plasticity state variable evolution rate component userobject class.
//
#include "CrystalPlasticityStateVariableEvolutionRateComponent.h"

template<>
InputParameters validParams<CrystalPlasticityStateVariableEvolutionRateComponent>()
{
  InputParameters params = validParams<CrystalPlasticityUOBase>();
  params.addRequiredParam<int >("nss", "Number of slip systems");
  params.addClassDescription("Crystal plasticity state variable evolution rate component base class.  Override the virtual functions in your class");
  return params;
}

CrystalPlasticityStateVariableEvolutionRateComponent::CrystalPlasticityStateVariableEvolutionRateComponent(const InputParameters & parameters) :
  CrystalPlasticityUOBase(parameters),
  _nss(getParam<int>("nss"))
{  
}

std::string
CrystalPlasticityStateVariableEvolutionRateComponent::crystalPlasticityUOName() const
{
  return "uo_state_evol_rate_comp_" + _uo_name;
}

