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
  params.addRequiredParam<Real>("apb_shear_energy", "Anti-phase boundary shear resistance energy");
  params.addRequiredParam<Real>("burgers_length", "Length of Burgers vector");
  params.addParam<Real>("number_density", 0.0, "Average number density of precipitate");
  params.addParam<Real>("size", 0.0, "Average size of precipitate");
  params.addParam<Real>("shear_mod", 79731, "Shear modulus in MPa");
  params.addClassDescription("Dislocation based constitutive mode userobject class for thermal slip resistance.  Override the virtual functions in your class");
  return params;
}

CPDislocationBasedAPBSlipResistance::CPDislocationBasedAPBSlipResistance(const InputParameters & parameters) :
    CrystalPlasticitySlipResistance(parameters),
    _apb_shear_energy(getParam<Real>("apb_shear_energy")),
    _b(getParam<Real>("burgers_length")),
    _number_density(getParam<Real>("number_density")),
    _size(getParam<Real>("size")),
    _shear_mod(getParam<Real>("shear_mod"))
{
}

bool
CPDislocationBasedAPBSlipResistance::calcSlipResistance(unsigned int qp, std::vector<Real> & val) const
{
  Real f = _number_density * 4.0/3.0 * libMesh::pi * std::pow(_size/2.0, 3.0);
  Real apb_shear_resist = 2.0 * _apb_shear_energy/ std::pow(_b, 2.0) * std::sqrt(f * _size/(libMesh::pi * _shear_mod));
  
  for (unsigned int i = 0; i < _variable_size; ++i)
    val[i] = apb_shear_resist;

  return true;
}
