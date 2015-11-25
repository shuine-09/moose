/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
//  Crystal plasticity slip resistance userobject class.
//
#include "CrystalPlasticitySlipResistance.h"

template<>
InputParameters validParams<CrystalPlasticitySlipResistance>()
{
  InputParameters params = validParams<CrystalPlasticityUOBase>();
  params.addRequiredParam<int >("nss", "Number of slip systems");
  params.addClassDescription("Crystal plasticity slip resistance base class.  Override the virtual functions in your class");
  return params;
}

CrystalPlasticitySlipResistance::CrystalPlasticitySlipResistance(const InputParameters & parameters) :
  CrystalPlasticityUOBase(parameters),
  _nss(getParam<int>("nss"))
{
}

void
CrystalPlasticitySlipResistance::initSlipSysProps(std::vector<Real> & val) const
{
}

std::vector<Real> 
CrystalPlasticitySlipResistance::getHardnessParams() const
{
  std::vector<Real> val;
  return val;
}

std::string
CrystalPlasticitySlipResistance::crystalPlasticityUOName() const
{
  return "uo_slip_resistance_" + _uo_name;
}
