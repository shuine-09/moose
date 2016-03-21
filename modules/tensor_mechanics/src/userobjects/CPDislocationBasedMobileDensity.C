/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "CPDislocationBasedMobileDensity.h"

template<>
InputParameters validParams<CPDislocationBasedMobileDensity>()
{
  InputParameters params = validParams<CrystalPlasticityStateVariable>();
  params.addParam<Real>("initial_value", "Initial values of mobile dislocation density");
  params.addClassDescription("Dislocation based constitutive model userobject class for mobile dislocation density.  Override the virtual functions in your class");
  return params;
}

CPDislocationBasedMobileDensity::CPDislocationBasedMobileDensity(const InputParameters & parameters) :
    CrystalPlasticityStateVariable(parameters),
    _initial_value(getParam<Real>("initial_value"))
{
}

void
CPDislocationBasedMobileDensity::initSlipSysProps(std::vector<Real> & val) const
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
CPDislocationBasedMobileDensity::assignSlipSysRes(std::vector<Real> & val) const
{
}

void
CPDislocationBasedMobileDensity::getInitSlipSysRes(std::vector<Real> & val) const
{
  if (_initial_value <= 0.0)
    mooseError("CPDislocationBasedMobileDensity: Initial value of mobile dislocation density is non positive");
  else
    for (unsigned int i = 0; i < _variable_size; ++i)
      val[i] = _initial_value;
}
