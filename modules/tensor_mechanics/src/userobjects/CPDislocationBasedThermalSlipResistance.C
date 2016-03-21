/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "CPDislocationBasedThermalSlipResistance.h"

template<>
InputParameters validParams<CPDislocationBasedThermalSlipResistance>()
{
  InputParameters params = validParams<CrystalPlasticitySlipResistance>();
  params.addRequiredParam<Real>("thermal_resist", "Thermal resistance to slip");
  params.addClassDescription("Dislocation based constitutive mode userobject class for thermal slip resistance.  Override the virtual functions in your class");
  return params;
}

CPDislocationBasedThermalSlipResistance::CPDislocationBasedThermalSlipResistance(const InputParameters & parameters) :
    CrystalPlasticitySlipResistance(parameters),
    _thermal_resist(getParam<Real>("thermal_resist"))
{
}

bool
CPDislocationBasedThermalSlipResistance::calcSlipResistance(unsigned int qp, std::vector<Real> & val) const
{
  for (unsigned int i = 0; i < _variable_size; ++i)
    val[i] = _thermal_resist;

  return true;
}
