/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "CPDislocationBasedMobileGlideRateComp.h"

template<>
InputParameters validParams<CPDislocationBasedMobileGlideRateComp>()
{
  InputParameters params = validParams<CrystalPlasticityStateVarRateComponent>();
  params.addParam<std::string>("uo_mobile_dislocation_density_name", "Name of mobile dislocation density: Same as state variable user object specified in input file.");
  params.addParam<std::string>("uo_immobile_dislocation_density_name", "Name of immobile dislocation density: Same as state variable user object specified in input file.");
  params.addParam<std::string>("uo_glide_slip_rate_name", "Name of glide slip rate: Same as state variable user object specified in input file.");
  params.addRequiredParam<Real>("burgers_length", "Length of Burgers vector");
  params.addParam<Real>("rho_mult_factor", 0.143, "Dislocation multiplication factor");
  params.addParam<Real>("rho_m_capture_radius", 1e-6, "Capture radius for mutual annihilation of mobile dislocations: Unit in mm");
  params.addParam<Real>("rho_imm_factor", 0.36, "Immobilization factor due to dislocations");
  params.addClassDescription("Dislocation based constitutive model userobject class for mobile dislocation density rate component for glide mechanism. Override the virtual functions in your class");
  return params;
}

CPDislocationBasedMobileGlideRateComp::CPDislocationBasedMobileGlideRateComp(const InputParameters & parameters) :
    CrystalPlasticityStateVarRateComponent(parameters),
    _mat_prop_mobile_dislocation_density(getMaterialProperty<std::vector<Real> >(parameters.get<std::string>("uo_mobile_dislocation_density_name"))),
    _mat_prop_immobile_dislocation_density(getMaterialProperty<std::vector<Real> >(parameters.get<std::string>("uo_immobile_dislocation_density_name"))),
    _mat_prop_glide_slip_rate(getMaterialProperty<std::vector<Real> >(parameters.get<std::string>("uo_glide_slip_rate_name"))),
    _b(getParam<Real>("burgers_length")),
    _k_mul(getParam<Real>("rho_mult_factor")),
    _r_c(getParam<Real>("rho_m_capture_radius")),
    _beta_rho(getParam<Real>("rho_imm_factor"))
{
}

bool
CPDislocationBasedMobileGlideRateComp::calcStateVariableEvolutionRateComponent(unsigned int qp, std::vector<Real> & val) const
{
  val.assign(_variable_size, 0.0);

  Real tot_rho_m = 0.0;
  for (unsigned int i = 0; i < _variable_size; ++i)
    tot_rho_m += _mat_prop_mobile_dislocation_density[qp][i];

  Real d_self, d_ann, d_imm, lambda_inv;

  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    d_self = _k_mul * std::sqrt(tot_rho_m) * std::abs(_mat_prop_glide_slip_rate[qp][i])/_b;
    d_ann = 2.0 * _r_c * _mat_prop_mobile_dislocation_density[qp][i] * std::abs(_mat_prop_glide_slip_rate[qp][i]) / _b;

    lambda_inv = _beta_rho * std::sqrt(_mat_prop_mobile_dislocation_density[qp][i] + _mat_prop_immobile_dislocation_density[qp][i]);
    d_imm = std::abs(_mat_prop_glide_slip_rate[qp][i]) * lambda_inv / _b ;

    val[i] = d_self - d_ann - d_imm;
  }

  return true;
}
