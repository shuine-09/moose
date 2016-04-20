/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "CPDislocationBasedGlideSlipRate.h"

template<>
InputParameters validParams<CPDislocationBasedGlideSlipRate>()
{
  InputParameters params = validParams<CrystalPlasticitySlipRate>();
  params.addParam<std::string>("uo_mobile_dislocation_density_name", "Name of mobile dislocation density: Same as state variable user object specified in input file.");
  params.addParam<std::string>("uo_thermal_slip_resistance_name", "Name of slip resistance property: Same as slip resistance user object specified in input file.");
  params.addParam<std::string>("uo_athermal_slip_resistance_name", "Name of slip resistance property: Same as slip resistance user object specified in input file.");
  params.addParam<Real>("penalty_param", 1e-3, "Penalty parameter for linear regularization");
  params.addRequiredParam<Real>("burgers_length", "Length of Burgers vector");
  params.addRequiredParam<Real>("disloc_line_dist", "Dislocation line distance");
  params.addRequiredParam<Real>("jump_freq", "Jump frequency");
  params.addRequiredParam<Real>("active_enthal", "Activation enthalpy");
  params.addParam<Real>("exponentp", 0.28, "Exponent p used in flow rule");
  params.addParam<Real>("exponentq", 1.34, "Exponent q used in flow rule");
  params.addParam<Real>("boltz_const", 1.38065e-20, "Boltzman constant: Default unit MPa-mm^3");
  params.addParam<Real>("temp", 273.0, "Temperature in K");
  params.addClassDescription("Dislocation based constitutive mode userobject class for glide slip rate in the matrix.  Override the virtual functions in your class");
  return params;
}

CPDislocationBasedGlideSlipRate::CPDislocationBasedGlideSlipRate(const InputParameters & parameters) :
    CrystalPlasticitySlipRate(parameters),
    _mat_prop_mobile_dislocation_density(getMaterialProperty<std::vector<Real> >(parameters.get<std::string>("uo_mobile_dislocation_density_name"))),
    _mat_prop_thermal_slip_resistance(getMaterialProperty<std::vector<Real> >(parameters.get<std::string>("uo_thermal_slip_resistance_name"))),
    _mat_prop_athermal_slip_resistance(getMaterialProperty<std::vector<Real> >(parameters.get<std::string>("uo_athermal_slip_resistance_name"))),
    _pk2(getMaterialPropertyByName<RankTwoTensor>("pk2")),
    _flow_direction(getMaterialProperty<std::vector<RankTwoTensor> >(_name + "_flow_direction")),
    _penalty_param(getParam<Real>("penalty_param")),
    _b(getParam<Real>("burgers_length")),
    _lg(getParam<Real>("disloc_line_dist")),
    _jump_freq(getParam<Real>("jump_freq")),
    _enthal(getParam<Real>("active_enthal")),
    _p(getParam<Real>("exponentp")),
    _q(getParam<Real>("exponentq")),
    _k(getParam<Real>("boltz_const")),
    _temp(getParam<Real>("temp"))
{
}

void
CPDislocationBasedGlideSlipRate::calcFlowDirection(unsigned int qp, std::vector<RankTwoTensor> & flow_direction) const
{
  DenseVector<Real> mo(LIBMESH_DIM*_variable_size),no(LIBMESH_DIM*_variable_size);

  // Update slip direction and normal with crystal orientation
  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
    {
      mo(i*LIBMESH_DIM+j) = 0.0;
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
        mo(i*LIBMESH_DIM+j) = mo(i*LIBMESH_DIM+j) + _crysrot[qp](j,k) * _mo(i*LIBMESH_DIM+k);
    }

    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
    {
      no(i*LIBMESH_DIM+j) = 0.0;
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
        no(i*LIBMESH_DIM+j) = no(i*LIBMESH_DIM+j) + _crysrot[qp](j,k) * _no(i*LIBMESH_DIM+k);
    }
  }

  // Calculate Schmid tensor and resolved shear stresses
  for (unsigned int i = 0; i < _variable_size; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
        flow_direction[i](j,k) = mo(i*LIBMESH_DIM+j) * no(i*LIBMESH_DIM+k);
}

bool
CPDislocationBasedGlideSlipRate::calcSlipRate(unsigned int qp, Real dt, std::vector<Real> & val) const
{
  DenseVector<Real> tau(_variable_size);

  for (unsigned int i = 0; i < _variable_size; ++i)
    tau(i) = _pk2[qp].doubleContraction(_flow_direction[qp][i]);

  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    Real sgn_tau = 1.0;
    if (tau(i) < 0.0)
      sgn_tau = -1.0;

    Real c = (std::abs(tau(i)) - _mat_prop_athermal_slip_resistance[qp][i])/_mat_prop_thermal_slip_resistance[qp][i];

    Real v1 = (std::abs(c) + c)/2.0;
    Real v2 = (std::abs(c) - c)/2.0;

    Real v1p = std::pow(v1, _p);
    if (v1p > 1.0)
    {
#ifdef DEBUG
      mooseWarning("CPDislocationBasedGlideSlipRate: Flow rule upper limit exceeded " << v1p);
#endif
      return false;
    }

    Real c2 = std::pow(1.0 - v1p, _q) + v2/_penalty_param;

    Real a = _mat_prop_mobile_dislocation_density[qp][i] * _b * _lg * _jump_freq;
    Real b = _enthal/ (_k * _temp);

    val[i] = a * std::exp(-b * c2) * sgn_tau;
  }

  return true;
}

bool
CPDislocationBasedGlideSlipRate::calcSlipRateDerivative(unsigned int qp, Real dt, std::vector<Real> & val) const
{
  DenseVector<Real> tau(_variable_size);

  for (unsigned int i = 0; i < _variable_size; ++i)
    tau(i) = _pk2[qp].doubleContraction(_flow_direction[qp][i]);

  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    Real sgn_tau = 1.0;
    if (tau(i) < 0.0)
      sgn_tau = -1.0;

    Real c = (std::abs(tau(i)) - _mat_prop_athermal_slip_resistance[qp][i])/_mat_prop_thermal_slip_resistance[qp][i];
    Real dc_dtau = sgn_tau / _mat_prop_thermal_slip_resistance[qp][i];

    Real sgn_c = 1.0;
    if (c < 0.0)
      sgn_c = -1.0;

    Real v1 = (std::abs(c) + c)/2.0;
    Real dv1_dc = (sgn_c + 1.0)/2.0;

    Real v2 = (std::abs(c) - c)/2.0;
    Real dv2_dc = (sgn_c - 1.0)/2.0;

    Real v1p = std::pow(v1, _p);

    Real c2 = std::pow(1.0 - v1p, _q) + v2/_penalty_param;

    Real dc2_dv1 = 0.0;

    if (v1 > 0.0)
      dc2_dv1 = - _q * _p * std::pow(1.0 - std::pow(v1, _p), _q - 1.0) * std::pow(v1, _p - 1.0);

    Real dc2_dv2 = 1.0/_penalty_param;

    Real dc2_dc = dc2_dv1 * dv1_dc + dc2_dv2 * dv2_dc;

    Real a = _mat_prop_mobile_dislocation_density[qp][i] * _b * _lg * _jump_freq;
    Real b = _enthal/ (_k * _temp);

    val[i] = -a * b * std::exp(-b * c2) * sgn_tau * dc2_dc * dc_dtau;
  }

  return true;
}
