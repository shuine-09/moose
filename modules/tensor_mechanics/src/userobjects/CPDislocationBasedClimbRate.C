/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "CPDislocationBasedClimbRate.h"

registerMooseObject("TensorMechanicsApp", CPDislocationBasedClimbRate);

template <>
InputParameters
validParams<CPDislocationBasedClimbRate>()
{
  InputParameters params = validParams<CrystalPlasticitySlipRate>();
  // params.addParam<std::string>("uo_concentration_name", "Name of concentration variable");
  params.addParam<std::string>("uo_mobile_dislocation_density_name",
                               "Name of mobile dislocation density: Same as "
                               "state variable user object specified in input "
                               "file.");
  params.addParam<std::string>("uo_immobile_dislocation_density_name",
                               "Name of immobile dislocation density: Same as "
                               "state variable user object specified in input "
                               "file.");
  params.addRequiredParam<Real>("diffusivity", "Self diffusivity");
  params.addRequiredParam<Real>("core_radius", "Dislocation core radius");
  params.addRequiredParam<Real>("burgers_length", "Length of Burgers vector");
  params.addRequiredParam<Real>("molar_volume", "Molar volume");
  params.addRequiredParam<Real>("atom_volume", "Atom volume");
  params.addRequiredParam<Real>("lattice_constant",
                                "Lattice constatnt: average inter-atomic distance.");
  params.addParam<Real>("boltz_const", "Boltzmann's constant");
  params.addParam<Real>("gas_constant", 8.314, "Gas constant J/mole");
  params.addParam<Real>("temp", 273.0, "Temperature in K");
  params.addParam<Real>("stress_factor", 1.0, "Stress concentration factor due to precipitate.");
  params.addParam<Real>("diffusivity_factor", 1.0, "diffusivity factor.");
  params.addRequiredParam<Real>("activation_energy", "Activation energy J/mole");
  params.addParam<Real>("precipitate_radius", 0.0, "The radius of precipitate size.");
  params.addRequiredParam<Real>("precipitate_volume_fraction",
                                "The precipitate_volume_fraction of precipitates.");
  params.addParam<Real>(
      "theta", libMesh::pi / 2, "The angle between dislocation line direction and burger vector.");

  params.addParam<bool>("use_damage", false, "Use damage.");
  params.addClassDescription("Dislocation based constitutive mode userobject "
                             "class for climb rate.  Override the virtual "
                             "functions in your class");
  return params;
}

CPDislocationBasedClimbRate::CPDislocationBasedClimbRate(const InputParameters & parameters)
  : CrystalPlasticitySlipRate(parameters),
    // _cv(getMaterialProperty<std::vector<Real>>(
    //     parameters.get<std::string>("uo_concentration_name"))),
    _rho_m(getMaterialProperty<std::vector<Real>>(
        parameters.get<std::string>("uo_mobile_dislocation_density_name"))),
    _rho_i(getMaterialProperty<std::vector<Real>>(
        parameters.get<std::string>("uo_immobile_dislocation_density_name"))),
    _pk2(getMaterialPropertyByName<RankTwoTensor>("pk2")),
    _flow_direction(getMaterialProperty<std::vector<RankTwoTensor>>(_name + "_flow_direction")),
    _diffusivity(getParam<Real>("diffusivity")),
    _rc(getParam<Real>("core_radius")),
    _b(getParam<Real>("burgers_length")),
    _molar_volume(getParam<Real>("molar_volume")),
    _atom_volume(getParam<Real>("atom_volume")),
    _lattice_constant(getParam<Real>("lattice_constant")),
    _boltz_const(getParam<Real>("boltz_const")),
    _gas_constant(getParam<Real>("gas_constant")),
    _temp(getParam<Real>("temp")),
    _stress_factor(getParam<Real>("stress_factor")),
    _diffusivity_factor(getParam<Real>("diffusivity_factor")),
    _activation_energy(getParam<Real>("activation_energy")),
    _precipitate_radius(getParam<Real>("precipitate_radius")),
    _precipitate_volume_fraction(getParam<Real>("precipitate_volume_fraction")),
    _theta(getParam<Real>("theta")),
    _use_damage(getParam<bool>("use_damage")),
    _w_old(_use_damage ? &getMaterialPropertyOldByName<Real>("w") : nullptr)
{
}

void
CPDislocationBasedClimbRate::calcFlowDirection(unsigned int qp,
                                               std::vector<RankTwoTensor> & flow_direction) const
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
CPDislocationBasedClimbRate::calcSlipRate(unsigned int qp, Real dt, std::vector<Real> & val) const
{
  DenseVector<Real> tau(_variable_size);

  const Real cv = 1.0;

  Real w = 0.0;

  if (_w_old)
    w = (*_w_old)[qp];

  for (unsigned int i = 0; i < _variable_size; ++i)
    tau(i) = -_pk2[qp].doubleContraction(_flow_direction[qp][i]) * _b / (1 - w);

  Real rho = 0.0;
  for (unsigned int i = 0; i < _variable_size; ++i)
    rho += (_rho_m[qp][i] + _rho_i[qp][i]);

  Real R = 1.0 / (2.0 * std::sqrt(rho));
  Real Dv_sd =
      _diffusivity * _diffusivity_factor * std::exp(-_activation_energy / (_gas_constant * _temp));
  // Real pre_factor = 12.0 * libMesh::pi * _atom_volume * Dv_sd /
  //                  (_b * (libMesh::pi * pow(_lattice_constant, 3.0) -
  //                         6.0 * _atom_volume * std::log(_lattice_constant / R)));
  Real pre_factor = 2.0 * libMesh::pi * Dv_sd / std::log(R / 4 / _b) / _b;

  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    // std::cout << "value = "
    //          << _stress_factor * tau(i) * _atom_volume /
    //                 (_boltz_const * _b * _temp)
    //<< std::endl;
    Real vc = -pre_factor / std::sin(_theta) *
              (std::exp(-_stress_factor * tau(i) * _atom_volume / (_boltz_const * _b * _temp) /
                        std::sin(_theta)) -
               1.0 * cv);
    // Real vc =
    //     pre_factor *
    //     (std::exp(-_stress_factor * tau(i) * _atom_volume / (_boltz_const * _b * _temp)) - 1.0);

    if (abs(vc) > 1.0e10)
    {
      vc = 0.0;
      // std::cout << "Dv_sd = " << Dv_sd << ", prefactor = "
      //           << std::exp(-_stress_factor * tau(i) * _atom_volume / (_boltz_const * _b *
      //           _temp))
      //          << ", tau = " << tau(i) << std::endl;
    }

    // val[i] = -_prefactor * (_rho_m[qp][i] + _rho_i[qp][i]) * _b * vc;
    val[i] = -_prefactor * _precipitate_volume_fraction * (_rho_m[qp][i]) * _b * vc;

    // std::cout << "slip system(" << i << ") : climb velocity = " << vc << ", climb rate = " <<
    // val[i]
    //           << ", resolved stress = " << tau(i) << " Pk2 yy = " << _pk2[qp](1, 1) << std::endl;

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
CPDislocationBasedClimbRate::calcSlipRateDerivative(unsigned int qp,
                                                    Real dt,
                                                    std::vector<Real> & val) const
{
  DenseVector<Real> tau(_variable_size);

  Real w = 0.0;
  if (_w_old)
    w = (*_w_old)[qp];

  for (unsigned int i = 0; i < _variable_size; ++i)
    tau(i) = -_pk2[qp].doubleContraction(_flow_direction[qp][i]) * _b / (1 - w);

  Real rho = 0.0;
  for (unsigned int i = 0; i < _variable_size; ++i)
    rho += (_rho_m[qp][i] + _rho_i[qp][i]);

  Real R = 1.0 / (2.0 * std::sqrt(rho));
  Real Dv_sd =
      _diffusivity * _diffusivity_factor * std::exp(-_activation_energy / (_gas_constant * _temp));
  // Real pre_factor = 12.0 * libMesh::pi * _atom_volume * Dv_sd /
  //                  (_b * (libMesh::pi * pow(_lattice_constant, 3.0) -
  //                         6.0 * _atom_volume * std::log(_lattice_constant / R)));
  Real pre_factor = 2.0 * libMesh::pi * Dv_sd / std::log(R / 4 / _b) / _b;

  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    Real dvc_dtau =
        pre_factor / std::sin(_theta) *
        (_stress_factor * _atom_volume / (_boltz_const * _b * _temp / std::sin(_theta)) *
         (std::exp(-_stress_factor * tau(i) * _atom_volume / (_boltz_const * _b * _temp))));

    // val[i] = -_prefactor * (_rho_m[qp][i] + _rho_i[qp][i]) * _b * dvc_dtau *
    //         (-1.0 * _b); // tau = -1.0 * pk2 * direction

    if (abs(dvc_dtau) > 1.0e10)
    {
      dvc_dtau = 0.0;
      // std::cout << "Dv_sd = " << Dv_sd << ", prefactor = "
      //           << std::exp(-_stress_factor * tau(i) * _atom_volume / (_boltz_const * _b *
      //           _temp))
      //           << ", tau = " << tau(i) << std::endl;
    }

    // val[i] = _prefactor * (_rho_i[qp][i] + _rho_m[qp][i]) * _b * dvc_dtau *
    //         (-1.0 * _b); // tau = -1.0 * pk2 * direction
    val[i] = -_prefactor * _precipitate_volume_fraction * (_rho_m[qp][i]) * _b * dvc_dtau *
             (-1.0 * _b); // tau = -1.0 * pk2 * direction
  }

  return true;
}
