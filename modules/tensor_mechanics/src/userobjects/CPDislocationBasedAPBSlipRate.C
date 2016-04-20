/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "CPDislocationBasedAPBSlipRate.h"

template<>
InputParameters validParams<CPDislocationBasedAPBSlipRate>()
{
  InputParameters params = validParams<CrystalPlasticitySlipRate>();
  params.addParam<std::string>("uo_mobile_dislocation_density_name", "Name of mobile dislocation density: Same as state variable user object specified in input file.");
  params.addParam<std::string>("uo_immobile_dislocation_density_name", "Name of immobile dislocation density: Same as state variable user object specified in input file.");
  params.addParam<std::string>("uo_thermal_slip_resistance_name", "Name of slip resistance property: Same as slip resistance user object specified in input file.");
  params.addParam<std::string>("uo_apb_slip_resistance_name", "Name of slip resistance property: Same as slip resistance user object specified in input file.");
  params.addRequiredParam<Real>("burgers_length", "Length of Burgers vector");
  params.addRequiredParam<Real>("disloc_line_dist", "Dislocation line distance");
  params.addRequiredParam<Real>("jump_freq", "Jump frequency");
  params.addRequiredParam<Real>("active_enthal", "Activation enthalpy");
  params.addParam<Real>("exponentp", 0.28, "Exponent p used in flow rule");
  params.addParam<Real>("exponentq", 1.34, "Exponent q used in flow rule");
  params.addParam<Real>("boltz_const", 1.38065e-20, "Boltzman constant: Default unit MPa-mm^3");
  params.addParam<Real>("temp", 273.0, "Temperature in K");
  params.addParam<Real>("factor", 1.0, "Scaling factor to capture the piling effect of dislocations at precipitate interface");
  params.addParam<Real>("critial_density", 0.0, "Critial dislocation density");
  params.addClassDescription("Dislocation based constitutive mode userobject class for APB shear slip rate in the matrix.  Override the virtual functions in your class");
  return params;
}

CPDislocationBasedAPBSlipRate::CPDislocationBasedAPBSlipRate(const InputParameters & parameters) :
    CrystalPlasticitySlipRate(parameters),
    _mat_prop_mobile_dislocation_density(getMaterialProperty<std::vector<Real> >(parameters.get<std::string>("uo_mobile_dislocation_density_name"))),
    _mat_prop_immobile_dislocation_density(getMaterialProperty<std::vector<Real> >(parameters.get<std::string>("uo_immobile_dislocation_density_name"))),
    _mat_prop_thermal_slip_resistance(getMaterialProperty<std::vector<Real> >(parameters.get<std::string>("uo_thermal_slip_resistance_name"))),
    _mat_prop_apb_slip_resistance(getMaterialProperty<std::vector<Real> >(parameters.get<std::string>("uo_apb_slip_resistance_name"))),
    _pk2(getMaterialPropertyByName<RankTwoTensor>("pk2")),
    _flow_direction(getMaterialProperty<std::vector<RankTwoTensor> >(_name + "_flow_direction")),
    _b(getParam<Real>("burgers_length")),
    _lg(getParam<Real>("disloc_line_dist")),
    _jump_freq(getParam<Real>("jump_freq")),
    _enthal(getParam<Real>("active_enthal")),
    _k(getParam<Real>("boltz_const")),
    _temp(getParam<Real>("temp")),
    _factor(getParam<Real>("factor")),
    _critial_density(getParam<Real>("critial_density"))
{
}

void
CPDislocationBasedAPBSlipRate::calcFlowDirection(unsigned int qp, std::vector<RankTwoTensor> & flow_direction) const
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
CPDislocationBasedAPBSlipRate::calcSlipRate(unsigned int qp, Real dt, std::vector<Real> & val) const
{
  DenseVector<Real> tau(_variable_size);

  for (unsigned int i = 0; i < _variable_size; ++i)
    tau(i) = _pk2[qp].doubleContraction(_flow_direction[qp][i]);

  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    Real c0 = _factor * _mat_prop_immobile_dislocation_density[qp][i] - _critial_density;
    if ( c0 <= 0.0 )
      continue;

    Real sgn_tau = 1.0;
    if (tau(i) < 0.0)
      sgn_tau = -1.0;

    Real c = (std::abs(tau(i)) -  _mat_prop_apb_slip_resistance[qp][i]);

    if ( c <= 0.0)
      c = 0.0;

    c = 1.0 - c / _mat_prop_thermal_slip_resistance[qp][i];

    Real a = _mat_prop_mobile_dislocation_density[qp][i] * _b * _lg * _jump_freq;
    Real b = _enthal/ (_k * _temp);

    val[i] = a * std::exp(-b * c) * sgn_tau * c0;
  }

  return true;
}

bool
CPDislocationBasedAPBSlipRate::calcSlipRateDerivative(unsigned int qp, Real dt, std::vector<Real> & val) const
{
  DenseVector<Real> tau(_variable_size);

  for (unsigned int i = 0; i < _variable_size; ++i)
    tau(i) = _pk2[qp].doubleContraction(_flow_direction[qp][i]);

  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    Real c0 = _factor * _mat_prop_immobile_dislocation_density[qp][i] - _critial_density;
    if ( c0 <= 0.0 )
      continue;

    Real sgn_tau = 1.0;
    if (tau(i) < 0.0)
      sgn_tau = -1.0;

    Real c = (std::abs(tau(i)) -  _mat_prop_apb_slip_resistance[qp][i]);
    Real dc_dtau = -sgn_tau / _mat_prop_thermal_slip_resistance[qp][i];

    if ( c <= 0.0)
      c = 0.0;

    c = 1.0 - c / _mat_prop_thermal_slip_resistance[qp][i];

    Real a = _mat_prop_mobile_dislocation_density[qp][i] * _b * _lg * _jump_freq;
    Real b = _enthal/ (_k * _temp);

    val[i] = -a * b * std::exp(-b * c) * sgn_tau * c0 * dc_dtau;
  }

  return true;
}
