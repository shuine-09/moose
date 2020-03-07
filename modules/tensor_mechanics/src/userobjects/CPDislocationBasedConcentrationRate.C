/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "CPDislocationBasedConcentrationRate.h"

registerMooseObject("TensorMechanicsApp", CPDislocationBasedConcentrationRate);

template <>
InputParameters
validParams<CPDislocationBasedConcentrationRate>()
{
  InputParameters params = validParams<CrystalPlasticityStateVarRateComponent>();
  params.addParam<std::string>("uo_concentration_name", "Name of concentration variable");
  params.addParam<std::string>("uo_climb_rate_name", "Name of climb rate variable");
  params.addParam<std::string>("uo_mobile_dislocation_density_name",
                               "Name of mobile dislocation density: Same as "
                               "state variable user object specified in input "
                               "file.");
  params.addParam<std::string>("uo_immobile_dislocation_density_name",
                               "Name of immobile dislocation density: Same as "
                               "state variable user object specified in input "
                               "file.");
  params.addRequiredParam<Real>("bulk_concentration", "Bulk concentration");
  params.addRequiredParam<Real>("diffusivity", "Self diffusivity");
  params.addRequiredParam<Real>("core_radius", "Dislocation core radius");
  params.addRequiredParam<Real>("burgers_length", "Length of Burgers vector");
  params.addRequiredParam<Real>("molar_volume", "Molar volume");
  params.addParam<Real>("gas_constant", 8.314, "Gas constant J/mole");
  params.addParam<Real>("temp", 273.0, "Temperature in K");
  params.addParam<Real>("boltz_const", "Boltzmann's constant");
  params.addParam<Real>("stress_factor", 1.0, "Stress concentration factor due to precipitate.");
  params.addParam<Real>("diffusivity_factor", 1.0, "diffusivity factor.");
  params.addRequiredParam<Real>("atom_volume", "Atom volume");
  params.addRequiredParam<Real>("activation_energy", "Activation energy J/mole");

  params.addClassDescription("Dislocation based constitutive mode userobject "
                             "class for climb rate.  Override the virtual "
                             "functions in your class");
  return params;
}

CPDislocationBasedConcentrationRate::CPDislocationBasedConcentrationRate(
    const InputParameters & parameters)
  : CrystalPlasticityStateVarRateComponent(parameters),
    _cv(getMaterialProperty<std::vector<Real>>(
        parameters.get<std::string>("uo_concentration_name"))),
    _flow_direction(getMaterialProperty<std::vector<RankTwoTensor>>(
        parameters.get<std::string>("uo_climb_rate_name") + "_flow_direction")),
    _rho_m(getMaterialProperty<std::vector<Real>>(
        parameters.get<std::string>("uo_mobile_dislocation_density_name"))),
    _rho_i(getMaterialProperty<std::vector<Real>>(
        parameters.get<std::string>("uo_immobile_dislocation_density_name"))),
    _pk2(getMaterialPropertyByName<RankTwoTensor>("pk2")),
    _c0(getParam<Real>("bulk_concentration")),
    _diffusivity(getParam<Real>("diffusivity")),
    _rc(getParam<Real>("core_radius")),
    _b(getParam<Real>("burgers_length")),
    _molar_volume(getParam<Real>("molar_volume")),
    _gas_constant(getParam<Real>("gas_constant")),
    _temp(getParam<Real>("temp")),
    _boltz_const(getParam<Real>("boltz_const")),
    _atom_volume(getParam<Real>("atom_volume")),
    _stress_factor(getParam<Real>("stress_factor")),
    _diffusivity_factor(getParam<Real>("diffusivity_factor")),
    _activation_energy(getParam<Real>("activation_energy"))
{
}

bool
CPDislocationBasedConcentrationRate::calcStateVariableEvolutionRateComponent(
    unsigned int qp, std::vector<Real> & val) const
{
  if (_variable_size > 1)
    mooseError("Concentration is a scalar variable of size 1");

  unsigned int size = _rho_m[qp].size();

  if (size != _rho_i[qp].size())
    mooseError("Dislocation density vectors should be of same length");

  if (size != _flow_direction[qp].size())
    mooseError("Flow direction and dislocation density vectors should be of "
               "same length");

  DenseVector<Real> tau(size);

  for (unsigned int i = 0; i < _variable_size; ++i)
    tau(i) = -_pk2[qp].doubleContraction(_flow_direction[qp][i]) * _b;

  Real rho = 0.0;
  for (unsigned int i = 0; i < size; ++i)
    rho += (_rho_m[qp][i] + _rho_i[qp][i]);

  Real R = 1.0 / (2.0 * std::sqrt(rho));
  Real Dv_sd =
      _diffusivity * _diffusivity_factor * std::exp(-_activation_energy / (_boltz_const * _temp));
  Real pre_factor = 2.0 * libMesh::pi * Dv_sd / std::log(R / 4 / _b) / _b;

  val[0] = 0.0;

  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    // std::cout << "value = "
    //          << _stress_factor * tau(i) * _atom_volume /
    //                 (_boltz_const * _b * _temp)
    //<< std::endl;
    Real vc = pre_factor *
              (std::exp(-_stress_factor * tau(i) * _atom_volume / (_boltz_const * _b * _temp)) -
               1.0 * _cv[qp][0]);

    // val[0] += -(_rho_m[qp][i] + _rho_i[qp][i]) * _b * vc;
    val[0] += -(_rho_i[qp][i]) * _b * vc;
  }

  // std::cout << "concentration rate = " << val[0] << std::endl;

  return true;
}
