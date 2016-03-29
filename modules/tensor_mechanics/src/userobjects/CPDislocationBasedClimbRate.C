/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "CPDislocationBasedClimbRate.h"

template<>
InputParameters validParams<CPDislocationBasedClimbRate>()
{
  InputParameters params = validParams<CrystalPlasticitySlipRate>();
  params.addParam<std::string>("uo_mobile_dislocation_density_name", "Name of mobile dislocation density: Same as state variable user object specified in input file.");
  params.addParam<std::string>("uo_immobile_dislocation_density_name", "Name of immobile dislocation density: Same as state variable user object specified in input file.");
  params.addRequiredParam<Real>("climb_resistance", "Climb resistance value");
  params.addRequiredParam<Real>("burgers_length", "Length of Burgers vector");
  params.addRequiredParam<Real>("mean_climb_distance", "Length of mean free path");
  params.addRequiredParam<Real>("exponent_p", "Exponent p");
  params.addRequiredParam<Real>("exponent_q", "Exponent q");
  params.addRequiredParam<Real>("active_enthal", "Activation enthalpy");
  params.addParam<Real>("boltz_const", 1.38065e-20, "Boltzman constant: Default unit MPa-mm^3");
  params.addParam<Real>("temp", 273.0, "Temperature in K");
  params.addParam<Real>("edge_screw_ratio", 1e8, "Ratio of edge to screw dislocation density");
  params.addClassDescription("Dislocation based constitutive mode userobject class for climb rate.  Override the virtual functions in your class");
  return params;
}

CPDislocationBasedClimbRate::CPDislocationBasedClimbRate(const InputParameters & parameters) :
    CrystalPlasticitySlipRate(parameters),
    _rho_m(getMaterialProperty<std::vector<Real> >(parameters.get<std::string>("uo_mobile_dislocation_density_name"))),
    _rho_i(getMaterialProperty<std::vector<Real> >(parameters.get<std::string>("uo_immobile_dislocation_density_name"))),
    _pk2(getMaterialPropertyByName<RankTwoTensor>("pk2")),
    _flow_direction(getMaterialProperty<std::vector<RankTwoTensor> >(_name + "_flow_direction")),
    _climb_resistance(getParam<Real>("climb_resistance")),
    _b(getParam<Real>("burgers_length")),
    _l(getParam<Real>("mean_climb_distance")),
    _p(getParam<Real>("exponent_p")),
    _q(getParam<Real>("exponent_q")),
    _enthal(getParam<Real>("active_enthal")),
    _k(getParam<Real>("boltz_const")),
    _temp(getParam<Real>("temp"))
{
  _edge_screw_angle = std::atan(getParam<Real>("edge_screw_ratio"));
}

void
CPDislocationBasedClimbRate::calcFlowDirection(unsigned int qp, std::vector<RankTwoTensor> & flow_direction) const
{
  RealVectorValue mo, no, to, ro;
  RankTwoTensor iden(RankTwoTensor::initIdentity);
  RankTwoTensor flow_dirn_tmp;

  // Update slip direction and normal with crystal orientation
  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
    {
      no(j) = 0.0;
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
        no(j) = no(j) + _crysrot[qp](j,k) * _no(i*LIBMESH_DIM+k);
    }
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
    {
      mo(j) = 0.0;
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
        mo(j) = mo(j) + _crysrot[qp](j,k) * _mo(i*LIBMESH_DIM+k);
    }
    
    to = no.cross(mo);
    ro = mo * std::sin(_edge_screw_angle) - to * std::cos(_edge_screw_angle);
    
    flow_dirn_tmp.vectorOuterProduct(mo, ro);
    flow_dirn_tmp = 0.5 * (flow_dirn_tmp + flow_dirn_tmp.transpose());
    flow_dirn_tmp = flow_dirn_tmp - flow_dirn_tmp.trace()/3.0 * iden;
    flow_direction[i] = flow_dirn_tmp;
  }
}

bool
CPDislocationBasedClimbRate::calcSlipRate(unsigned int qp, Real dt, std::vector<Real> & val) const
{
  DenseVector<Real> tau(_variable_size);

  for (unsigned int i = 0; i < _variable_size; ++i)
    tau(i) = _pk2[qp].doubleContraction(_flow_direction[qp][i]);

  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    Real tau_ratio = std::abs(tau(i)) / _climb_resistance;
    Real factor = std::pow(tau_ratio, _p);

    if (factor > 1.0)
    {
#ifdef DEBUG
      mooseWarning("Climb increment exceed range " << factor);
#endif
      return false;
    }

    val[i] = (_rho_m[qp][i] + _rho_i[qp][i]) * _b * _l * std::exp(-_enthal/(_k * _temp) * std::pow(1.0 - std::pow(tau_ratio, _p), _q)) * copysign(1.0, tau(i));

    if (std::abs(val[i] * dt) > _slip_incr_tol)
    {
#ifdef DEBUG
      mooseWarning("Maximum allowable climb increment exceeded " << std::abs(val[i])*dt);
#endif
      return false;
    }
  }
  return true;
}

bool
CPDislocationBasedClimbRate::calcSlipRateDerivative(unsigned int qp, Real dt, std::vector<Real> & val) const
{
  DenseVector<Real> tau(_variable_size);

  for (unsigned int i = 0; i < _variable_size; ++i)
    tau(i) = _pk2[qp].doubleContraction(_flow_direction[qp][i]);

  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    Real tau_ratio = std::abs(tau(i)) / _climb_resistance;
    Real factor = std::pow(tau_ratio, _p);

    if (factor > 1.0)
    {
#ifdef DEBUG
      mooseWarning("Climb increment exceed range " << factor);
#endif
      return false;
    }

    val[i] = (_rho_m[qp][i] + _rho_i[qp][i]) * _b * _l * std::exp(-_enthal/(_k * _temp) * std::pow(1.0 - factor, _q));

    if (std::abs(val[i] * dt) > _slip_incr_tol)
    {
#ifdef DEBUG
      mooseWarning("Maximum allowable climb increment exceeded " << std::abs(val[i])*dt);
#endif
      return false;
    }

    Real a1 = 1.0/_climb_resistance;
    Real a2 = _p * std::pow(tau_ratio, _p - 1.0);
    Real a3 = _q * std::pow(1.0 - factor, _q - 1.0);
    
    val[i] *= _enthal/(_k * _temp) * a1 * a2 * a3;
  }
  return true;
}
