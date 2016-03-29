/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "CPDislocationBasedAPBSlipResistance.h"

template<>
InputParameters validParams<CPDislocationBasedAPBSlipResistance>()
{
  InputParameters params = validParams<CrystalPlasticitySlipResistance>();
  params.addRequiredParam<Real>("apb_shear_resist", "Anti-phase boundary shear resistance to slip");
  params.addClassDescription("Dislocation based constitutive mode userobject class for thermal slip resistance.  Override the virtual functions in your class");
  return params;
}

CPDislocationBasedAPBSlipResistance::CPDislocationBasedAPBSlipResistance(const InputParameters & parameters) :
    CrystalPlasticitySlipResistance(parameters),
    _apb_shear_resist(getParam<Real>("apb_shear_resist"))
{
}

bool
CPDislocationBasedAPBSlipResistance::calcSlipResistance(unsigned int qp, std::vector<Real> & val) const
{
  for (unsigned int i = 0; i < _variable_size; ++i)
    val[i] = _apb_shear_resist;

  return true;
}
