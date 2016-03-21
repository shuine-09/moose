/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "CPDislocationBasedImmobileGlideRateComp.h"

template<>
InputParameters validParams<CPDislocationBasedImmobileGlideRateComp>()
{
  InputParameters params = validParams<CrystalPlasticityStateVarRateComponent>();
  params.addParam<std::string>("uo_mobile_dislocation_density_name", "Name of mobile dislocation density: Same as state variable user object specified in input file.");
  params.addParam<std::string>("uo_immobile_dislocation_density_name", "Name of immobile dislocation density: Same as state variable user object specified in input file.");
  params.addParam<std::string>("uo_glide_slip_rate_name", "Name of glide slip rate: Same as state variable user object specified in input file.");
  params.addRequiredParam<Real>("burgers_length", "Length of Burgers vector");
  params.addParam<Real>("rho_imm_factor", 0.36, "Immobilization factor due to dislocations");
  params.addClassDescription("Dislocation based constitutive model userobject class for mobile dislocation density rate component for glide mechanism. Override the virtual functions in your class");
  return params;
}

CPDislocationBasedImmobileGlideRateComp::CPDislocationBasedImmobileGlideRateComp(const InputParameters & parameters) :
    CrystalPlasticityStateVarRateComponent(parameters),
    _mat_prop_mobile_dislocation_density(getMaterialProperty<std::vector<Real> >(parameters.get<std::string>("uo_mobile_dislocation_density_name"))),
    _mat_prop_immobile_dislocation_density(getMaterialProperty<std::vector<Real> >(parameters.get<std::string>("uo_immobile_dislocation_density_name"))),
    _mat_prop_glide_slip_rate(getMaterialProperty<std::vector<Real> >(parameters.get<std::string>("uo_glide_slip_rate_name"))),
    _b(getParam<Real>("burgers_length")),
    _beta_rho(getParam<Real>("rho_imm_factor"))
{
}

bool
CPDislocationBasedImmobileGlideRateComp::calcStateVariableEvolutionRateComponent(unsigned int qp, std::vector<Real> & val) const
{
  val.assign(_variable_size, 0.0);

  Real lambda_inv = 0.0;

  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    lambda_inv = _beta_rho * std::sqrt(_mat_prop_mobile_dislocation_density[qp][i] + _mat_prop_immobile_dislocation_density[qp][i]);
    val[i] = std::abs(_mat_prop_glide_slip_rate[qp][i]) * lambda_inv / _b;
  }

  return true;
}
