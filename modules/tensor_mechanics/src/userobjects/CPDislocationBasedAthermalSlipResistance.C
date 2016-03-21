/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "CPDislocationBasedAthermalSlipResistance.h"

template<>
InputParameters validParams<CPDislocationBasedAthermalSlipResistance>()
{
  InputParameters params = validParams<CrystalPlasticitySlipResistance>();
  params.addParam<std::string>("uo_mobile_dislocation_density_name", "Name of mobile dislocation density: Same as state variable user object specified in input file.");
  params.addParam<std::string>("uo_immobile_dislocation_density_name", "Name of immobile dislocation density: Same as state variable user object specified in input file.");
  params.addRequiredParam<Real>("burgers_length", "Length of Burgers vector");
  params.addParam<Real>("rho_m_barrier_factor", 0.3, "Mobile dislocation barrier strength factor");
  params.addParam<Real>("shear_mod", 79731, "Shear modulus in MPa");
  params.addParam<Real>("rho_self_hard_factor", 1.0, "Self hardening of dislocations");
  params.addParam<Real>("rho_latent_hard_factor", 0.2, "Latent hardening factor");
  params.addClassDescription("Dislocation based constitutive mode userobject class for athermal slip resistance.  Override the virtual functions in your class");
  return params;
}

CPDislocationBasedAthermalSlipResistance::CPDislocationBasedAthermalSlipResistance(const InputParameters & parameters) :
    CrystalPlasticitySlipResistance(parameters),
    _mat_prop_mobile_dislocation_density(getMaterialProperty<std::vector<Real> >(parameters.get<std::string>("uo_mobile_dislocation_density_name"))),
    _mat_prop_immobile_dislocation_density(getMaterialProperty<std::vector<Real> >(parameters.get<std::string>("uo_immobile_dislocation_density_name"))),
    _b(getParam<Real>("burgers_length")),
    _q_p(getParam<Real>("rho_m_barrier_factor")),
    _shear_mod(getParam<Real>("shear_mod")),
    _self_harden(getParam<Real>("rho_self_hard_factor")),
    _latent_harden(getParam<Real>("rho_latent_hard_factor"))
{
  _interaction_matrix.resize( _variable_size * _variable_size );

  for (unsigned int i = 0; i < _variable_size; ++i)
    for (unsigned int j = 0; j < _variable_size; ++j)
    {
      _interaction_matrix[i * _variable_size + j] = _latent_harden;
      if (i == j) _interaction_matrix[i * _variable_size + j] = _self_harden;
    }
}

bool
CPDislocationBasedAthermalSlipResistance::calcSlipResistance(unsigned int qp, std::vector<Real> & val) const
{
  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    Real disloc = 0.0;
    for (unsigned int j = 0; j < _variable_size; ++j)
      disloc += _interaction_matrix[i * _variable_size + j] * (_mat_prop_mobile_dislocation_density[qp][j] + _mat_prop_immobile_dislocation_density[qp][j]);
    disloc *= _q_p;
    val[i] = _shear_mod * _b * std::sqrt(disloc);
  }

  return true;
}
