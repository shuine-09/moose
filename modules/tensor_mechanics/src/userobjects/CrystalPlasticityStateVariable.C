/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
//  Crystal plasticity state variable userobject class.
//
#include "CrystalPlasticityStateVariable.h"

template<>
InputParameters validParams<CrystalPlasticityStateVariable>()
{
  InputParameters params = validParams<CrystalPlasticityUOBase>();
  params.addRequiredParam<int >("nss", "Number of slip systems");
  params.addClassDescription("Crystal plasticity state variable class.  Override the virtual functions in your class");
  return params;
}

CrystalPlasticityStateVariable::CrystalPlasticityStateVariable(const InputParameters & parameters) :
  CrystalPlasticityUOBase(parameters),
  _nss(getParam<int>("nss"))
{}

std::string
CrystalPlasticityStateVariable::crystalPlasticityUOName() const
{
  return "uo_state_var_" + _uo_name;
}

