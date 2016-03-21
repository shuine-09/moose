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
  params.addParam<std::string>("uo_apb_slip_resistance_name", "Name of slip resistance property: Same as slip resistance user object specified in input file.");
  params.addRequiredParam<Real>("burgers_length", "Length of Burgers vector");
  params.addRequiredParam<Real>("disloc_line_dist", "Dislocation line distance");
  params.addRequiredParam<Real>("reference_shear_rate", "Reference shear rate");
  params.addParam<Real>("exponentm", 1.0, "Exponent m used in flow rule");
  params.addParam<FunctionName>("factor_function", "Function to obtain shear rate dependent factor");
  params.addClassDescription("Dislocation based constitutive mode userobject class for APB shear slip rate in the matrix.  Override the virtual functions in your class");
  return params;
}

CPDislocationBasedAPBSlipRate::CPDislocationBasedAPBSlipRate(const InputParameters & parameters) :
    CrystalPlasticitySlipRate(parameters),
    _mat_prop_mobile_dislocation_density(getMaterialProperty<std::vector<Real> >(parameters.get<std::string>("uo_mobile_dislocation_density_name"))),
    _mat_prop_immobile_dislocation_density(getMaterialProperty<std::vector<Real> >(parameters.get<std::string>("uo_immobile_dislocation_density_name"))),
    _mat_prop_apb_slip_resistance(getMaterialProperty<std::vector<Real> >(parameters.get<std::string>("uo_apb_slip_resistance_name"))),
    _pk2(getMaterialPropertyByName<RankTwoTensor>("pk2")),
    _flow_direction(getMaterialProperty<std::vector<RankTwoTensor> >(_name + "_flow_direction")),
    _b(getParam<Real>("burgers_length")),
    _lg(getParam<Real>("disloc_line_dist")),
    _a0(getParam<Real>("reference_shear_rate")),
    _m(getParam<Real>("exponentm")),
    _factor_function(isParamValid("factor_function") ? &getFunction("factor_function") : NULL)
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
  Point p;
  DenseVector<Real> tau(_variable_size);
  
  for (unsigned int i = 0; i < _variable_size; ++i)
    tau(i) = _pk2[qp].doubleContraction(_flow_direction[qp][i]);

  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    Real c = 0.0;
    Real scaled_tau = _factor_function->value(_mat_prop_immobile_dislocation_density[qp][i], p) * tau(i);

    if (_mat_prop_apb_slip_resistance[qp][i] > 0.0 && std::abs(scaled_tau) > _mat_prop_apb_slip_resistance[qp][i])
      c = std::pow(std::abs(scaled_tau)/_mat_prop_apb_slip_resistance[qp][i] - 1.0, _m);

    Real a = _mat_prop_mobile_dislocation_density[qp][i] * _b * _lg * _a0 * c;

    Real sgn_tau = 1.0;
    if (scaled_tau < 0.0)
      sgn_tau = -1.0;

    val[i] = _prefactor * a * sgn_tau;
  }

  return true;
}

bool
CPDislocationBasedAPBSlipRate::calcSlipRateDerivative(unsigned int qp, Real dt, std::vector<Real> & val) const
{
  DenseVector<Real> tau(_variable_size);
  Point p;

  for (unsigned int i = 0; i < _variable_size; ++i)
    tau(i) = _pk2[qp].doubleContraction(_flow_direction[qp][i]);

  for (unsigned int i = 0; i < _variable_size; ++i)
  {

    Real dc_dtau = 0.0;
    Real factor = _factor_function->value(_mat_prop_immobile_dislocation_density[qp][i], p);
    Real scaled_tau = factor * tau(i);

    if (_mat_prop_apb_slip_resistance[qp][i] > 0.0 && std::abs(scaled_tau) > _mat_prop_apb_slip_resistance[qp][i])
      dc_dtau = _m * std::pow(std::abs(scaled_tau)/_mat_prop_apb_slip_resistance[qp][i] - 1.0, _m - 1.0)/_mat_prop_apb_slip_resistance[qp][i] * factor;

    Real a = _mat_prop_mobile_dislocation_density[qp][i] * _b * _lg * _a0 * dc_dtau;
    val[i] = _prefactor * a;
  }

  return true;
}
