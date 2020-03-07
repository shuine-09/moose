/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "CPDislocationBasedClimbRateDyson.h"

registerMooseObject("TensorMechanicsApp", CPDislocationBasedClimbRateDyson);

template <>
InputParameters
validParams<CPDislocationBasedClimbRateDyson>()
{
  InputParameters params = validParams<CrystalPlasticitySlipRate>();
  params.addParam<std::string>("uo_concentration_name", "Name of concentration variable");
  params.addParam<std::string>("uo_mobile_dislocation_density_name",
                               "Name of mobile dislocation density: Same as "
                               "state variable user object specified in input "
                               "file.");
  params.addParam<std::string>("uo_immobile_dislocation_density_name",
                               "Name of immobile dislocation density: Same as "
                               "state variable user object specified in input "
                               "file.");
  params.addRequiredParam<Real>("diffusivity", "Self diffusivity");
  params.addRequiredParam<Real>("burgers_length", "Length of Burgers vector");
  params.addParam<Real>("boltz_const", "Boltzmann's constant");
  params.addParam<Real>("temp", 273.0, "Temperature in K");
  params.addParam<Real>("stress_factor", 1.0, "Stress concentration factor due to precipitate.");
  params.addParam<Real>("diffusivity_factor", 1.0, "diffusivity factor.");
  params.addRequiredParam<Real>("activation_energy", "Activation energy J/mole");
  params.addRequiredParam<Real>("precipitate_radius", "The radius of precipitate size.");
  params.addRequiredParam<Real>("precipitate_volume_fraction",
                                "The precipitate_volume_fraction of precipitates.");
  params.addParam<Real>(
      "theta", libMesh::pi, "The angle between dislocation line direction and burger vector.");
  params.addClassDescription("Dislocation based constitutive mode userobject "
                             "class for climb rate.  Override the virtual "
                             "functions in your class");
  return params;
}

CPDislocationBasedClimbRateDyson::CPDislocationBasedClimbRateDyson(
    const InputParameters & parameters)
  : CrystalPlasticitySlipRate(parameters),
    _cv(getMaterialProperty<std::vector<Real>>(
        parameters.get<std::string>("uo_concentration_name"))),
    _rho_m(getMaterialProperty<std::vector<Real>>(
        parameters.get<std::string>("uo_mobile_dislocation_density_name"))),
    _rho_i(getMaterialProperty<std::vector<Real>>(
        parameters.get<std::string>("uo_immobile_dislocation_density_name"))),
    _pk2(getMaterialPropertyByName<RankTwoTensor>("pk2")),
    _flow_direction(getMaterialProperty<std::vector<RankTwoTensor>>(_name + "_flow_direction")),
    _diffusivity(getParam<Real>("diffusivity")),
    _b(getParam<Real>("burgers_length")),
    _boltz_const(getParam<Real>("boltz_const")),
    _temp(getParam<Real>("temp")),
    _stress_factor(getParam<Real>("stress_factor")),
    _diffusivity_factor(getParam<Real>("diffusivity_factor")),
    _activation_energy(getParam<Real>("activation_energy")),
    _precipitate_radius(getParam<Real>("precipitate_radius")),
    _precipitate_volume_fraction(getParam<Real>("precipitate_volume_fraction")),
    _theta(getParam<Real>("theta"))
{
}

void
CPDislocationBasedClimbRateDyson::calcFlowDirection(
    unsigned int qp, std::vector<RankTwoTensor> & flow_direction) const
{
  RealVectorValue mo;

  // Update slip direction and normal with crystal orientation
  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
    {
      mo(j) = 0.0;
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
        mo(j) = mo(j) + _crysrot[qp](j, k) * _mo(i * LIBMESH_DIM + k);
    }
    flow_direction[i].vectorOuterProduct(mo, mo);
  }
}

bool
CPDislocationBasedClimbRateDyson::calcSlipRate(unsigned int qp,
                                               Real dt,
                                               std::vector<Real> & val) const
{
  DenseVector<Real> tau(_variable_size);

  for (unsigned int i = 0; i < _variable_size; ++i)
    tau(i) = -_pk2[qp].doubleContraction(_flow_direction[qp][i]) * _b;

  Real Dv_sd =
      _diffusivity * _diffusivity_factor * std::exp(-_activation_energy / (_boltz_const * _temp));

  Real pre_factor = Dv_sd / _b;

  // Real spacing_precipitate =
  //     std::sqrt(8.0 / (libMesh::pi * 3.0 * _precipitate_volume_fraction)) * _precipitate_radius -
  //     _precipitate_radius;
  Real spacing_precipitate = _b;

  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    Real vc = -pre_factor * std::sinh(-_stress_factor * tau(i) * _b * spacing_precipitate /
                                      (_boltz_const * _temp));

    if (abs(vc) > 1.0e10)
      vc = 0.0;

    val[i] = -_prefactor * _precipitate_volume_fraction * (_rho_m[qp][i]) * _b * vc;

    if (std::abs(val[i] * dt) > _slip_incr_tol)
    {
#ifdef DEBUG
      mooseWarning("Maximum allowable climb increment exceeded ", std::abs(val[i]) * dt);
#endif
      return false;
    }
  }
  return true;
}

bool
CPDislocationBasedClimbRateDyson::calcSlipRateDerivative(unsigned int qp,
                                                         Real dt,
                                                         std::vector<Real> & val) const
{
  DenseVector<Real> tau(_variable_size);

  for (unsigned int i = 0; i < _variable_size; ++i)
    tau(i) = -_pk2[qp].doubleContraction(_flow_direction[qp][i]) * _b;

  Real Dv_sd =
      _diffusivity * _diffusivity_factor * std::exp(-_activation_energy / (_boltz_const * _temp));
  Real pre_factor = Dv_sd / _b;

  // Real spacing_precipitate =
  //     std::sqrt(8.0 / (libMesh::pi * 3.0 * _precipitate_volume_fraction)) * _precipitate_radius -
  //     _precipitate_radius;

  Real spacing_precipitate = _b;

  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    Real dvc_dtau =
        pre_factor *
        std::cosh(-_stress_factor * tau(i) * _b * spacing_precipitate / (_boltz_const * _temp)) *
        (_stress_factor * tau(i) * _b * spacing_precipitate / (_boltz_const * _temp));

    if (abs(dvc_dtau) > 1.0e10)
      dvc_dtau = 0.0;

    val[i] =
        -_prefactor * _precipitate_volume_fraction * (_rho_m[qp][i]) * _b * dvc_dtau * (-1.0 * _b);
  }

  return true;
}
