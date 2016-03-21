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
  params.addParam<std::string>("uo_concentration_name", "Name of concentration variable");
  params.addParam<std::string>("uo_mobile_dislocation_density_name", "Name of mobile dislocation density: Same as state variable user object specified in input file.");
  params.addParam<std::string>("uo_immobile_dislocation_density_name", "Name of immobile dislocation density: Same as state variable user object specified in input file.");
  params.addRequiredParam<Real>("bulk_concentration", "Bulk concentration");
  params.addRequiredParam<Real>("diffusivity", "Self diffusivity");
  params.addRequiredParam<Real>("core_radius", "Dislocation core radius");
  params.addRequiredParam<Real>("burgers_length", "Length of Burgers vector");
  params.addRequiredParam<Real>("molar_volume", "Molar volume");
  params.addParam<Real>("gas_constant", 8.314, "Gas constant J/mole");
  params.addParam<Real>("temp", 273.0, "Temperature in K");
  
  params.addClassDescription("Dislocation based constitutive mode userobject class for climb rate.  Override the virtual functions in your class");
  return params;
}

CPDislocationBasedClimbRate::CPDislocationBasedClimbRate(const InputParameters & parameters) :
    CrystalPlasticitySlipRate(parameters),
    _cv(getMaterialProperty<std::vector<Real> >(parameters.get<std::string>("uo_concentration_name"))),
    _rho_m(getMaterialProperty<std::vector<Real> >(parameters.get<std::string>("uo_mobile_dislocation_density_name"))),
    _rho_i(getMaterialProperty<std::vector<Real> >(parameters.get<std::string>("uo_immobile_dislocation_density_name"))),
    _pk2(getMaterialPropertyByName<RankTwoTensor>("pk2")),
    _flow_direction(getMaterialProperty<std::vector<RankTwoTensor> >(_name + "_flow_direction")),
    _c0(getParam<Real>("bulk_concentration")),
    _diffusivity(getParam<Real>("diffusivity")),
    _rc(getParam<Real>("core_radius")),
    _b(getParam<Real>("burgers_length")),
    _molar_volume(getParam<Real>("molar_volume")),
    _gas_constant(getParam<Real>("gas_constant")),
    _temp(getParam<Real>("temp"))
{
}

void
CPDislocationBasedClimbRate::calcFlowDirection(unsigned int qp, std::vector<RankTwoTensor> & flow_direction) const
{
  RealVectorValue mo;

  // Update slip direction and normal with crystal orientation
  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
    {
      mo(j) = 0.0;
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
        mo(j) = mo(j) + _crysrot[qp](j,k) * _mo(i*LIBMESH_DIM+k);
    }
    flow_direction[i].vectorOuterProduct(mo, mo);
  }
}

bool
CPDislocationBasedClimbRate::calcSlipRate(unsigned int qp, Real dt, std::vector<Real> & val) const
{
  DenseVector<Real> tau(_variable_size);

  for (unsigned int i = 0; i < _variable_size; ++i)
    tau(i) = -_pk2[qp].doubleContraction(_flow_direction[qp][i]);

  Real rho = 0.0;
  for (unsigned int i = 0; i < _variable_size; ++i)
    rho += (_rho_m[qp][i] + _rho_i[qp][i]);

  Real rinf = 1.0/(2.0 * std::sqrt(rho));
  Real pre_factor = 2.0 * libMesh::pi * _diffusivity/(_b * std::log(rinf/_rc));
  
  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    Real ceq = _c0 * std::exp(-tau(i) * _molar_volume/(_gas_constant * _temp));
    Real vc = pre_factor * (ceq - _cv[qp][0]);

    val[i] = _prefactor * (_rho_m[qp][i] + _rho_i[qp][i]) * _b * vc;

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

  Real rho = 0.0;
  for (unsigned int i = 0; i < _variable_size; ++i)
    rho += (_rho_m[qp][i] + _rho_i[qp][i]);

  Real rinf = 1.0/(2.0 * std::sqrt(rho));
  Real pre_factor = 2.0 * libMesh::pi * _diffusivity/(_b * std::log(rinf/_rc));
  
  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    Real dceq_dtau = -_c0 * std::exp(-tau(i) * _molar_volume/(_gas_constant * _temp)) * _molar_volume/(_gas_constant * _temp);
    Real dvc_dtau = pre_factor * dceq_dtau;

    val[i] = _prefactor * (_rho_m[qp][i] + _rho_i[qp][i]) * _b * dvc_dtau;
  }
  return true;
}
