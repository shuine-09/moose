/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "CPDislocationBasedImmobileDensity.h"

template<>
InputParameters validParams<CPDislocationBasedImmobileDensity>()
{
  InputParameters params = validParams<CrystalPlasticityStateVariable>();
  params.addParam<Real>("initial_value", "Initial values of immobile dislocation density");
  params.addClassDescription("Dislocation based constitutive model userobject class for immobile dislocation density.  Override the virtual functions in your class");
  return params;
}

CPDislocationBasedImmobileDensity::CPDislocationBasedImmobileDensity(const InputParameters & parameters) :
    CrystalPlasticityStateVariable(parameters),
    _initial_value(getParam<Real>("initial_value"))
{
}

void
CPDislocationBasedImmobileDensity::initSlipSysProps(std::vector<Real> & val) const
{
  switch (_intvar_read_type)
  {
    case 0:
      assignSlipSysRes(val);
      break;
    case 1:
      readFileInitSlipSysRes(val);
      break;
    default:
      getInitSlipSysRes(val);
  }
}

void
CPDislocationBasedImmobileDensity::assignSlipSysRes(std::vector<Real> & val) const
{
}

void
CPDislocationBasedImmobileDensity::getInitSlipSysRes(std::vector<Real> & val) const
{
  if (_initial_value <= 0.0)
    mooseError("CPDislocationBasedImmobileDensity: Initial value of immobile dislocation density is non positive");
  else
    for (unsigned int i = 0; i < _variable_size; ++i)
      val[i] = _initial_value;
}
