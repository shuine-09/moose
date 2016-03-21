/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "CPDislocationAndPptAthermalResistance.h"

template<>
InputParameters validParams<CPDislocationAndPptAthermalResistance>()
{
  InputParameters params = validParams<CPDislocationBasedAthermalSlipResistance>();
  params.addParam<std::string>("uo_shear_rate_name", "Name of shear rate property");
  params.addRequiredParam<Real>("apb_shear_energy", "Anti-phase boundary shear resistance energy");
  params.addParam<Real>("number_density", 0.0, "Average number density of precipitate");
  params.addParam<Real>("size", 0.0, "Average size of precipitate");
  params.addParam<FunctionName>("factor_function", "Function to obtain shear rate dependent factor");
  return params;
}

CPDislocationAndPptAthermalResistance::CPDislocationAndPptAthermalResistance(const InputParameters & parameters) :
    CPDislocationBasedAthermalSlipResistance(parameters),
    _shear_rate_prop(getMaterialProperty<std::vector<Real> >(parameters.get<std::string>("uo_shear_rate_name"))),
    _apb_shear_energy(getParam<Real>("apb_shear_energy")),
    _number_density(getParam<Real>("number_density")),
    _size(getParam<Real>("size")),
    _factor_function(isParamValid("factor_function") ? &getFunction("factor_function") : NULL)
{
}

bool
CPDislocationAndPptAthermalResistance::calcSlipResistance(unsigned int qp, std::vector<Real> & val) const
{
  Real f = _number_density * 4.0/3.0 * libMesh::pi * std::pow(_size/2.0, 3.0);
  Real apb_shear_resist = 2.0 * _apb_shear_energy/ std::pow(_b, 2.0) * std::sqrt(f * _size/(libMesh::pi * _shear_mod));

  Point p;
  if (CPDislocationBasedAthermalSlipResistance::calcSlipResistance(qp, val))
  {
    for (unsigned int i = 0; i < _variable_size; ++i)
    {
      Real disloc_resist = std::pow(val[i], 2.0);
      val[i] = disloc_resist + _factor_function->value(_shear_rate_prop[qp][i], p) * std::pow(apb_shear_resist, 2.0);
      val[i] = std::sqrt(val[i]); 
    }
  }
  else
    return false;

  return true;
}
