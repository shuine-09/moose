/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
//  Crystal plasticity state variable userobject class.
//
#include "CrystalPlasticityStateVariableGSS.h"

template<>
InputParameters validParams<CrystalPlasticityStateVariableGSS>()
{
  InputParameters params = validParams<CrystalPlasticityStateVariable>();
  params.addParam<bool>("use_slip_resistance_as_state_var", false, "Slip resistance is used as state varibale");
  params.addParam<UserObjectName>("slip_resistance_uo","The Slip resistance user object");
  params.addParam<std::vector<std::string> >("uo_state_var_evol_rate_name", "Name of state variable evolution rate property: Same as state variable evolution rate user object specified in input file.");
  params.addParam<std::vector<std::string> >("uo_state_var_name", "Name of state variable property: Same as state variable user object specified in input file.");
  params.addClassDescription("Crystal plasticity state class.  Override the virtual functions in your class");
  return params;
}

CrystalPlasticityStateVariableGSS::CrystalPlasticityStateVariableGSS(const InputParameters & parameters) :
  CrystalPlasticityStateVariable(parameters),
  _use_slip_resistance_as_state_var(getParam<bool>("use_slip_resistance_as_state_var")),
  _uo_slip_resistance(isParamValid("slip_resistance_uo") ? & getUserObject<CrystalPlasticitySlipResistance>("slip_resistance_uo") : NULL),
  _num_mat_state_var_evol_rates(parameters.get<std::vector<std::string> >("uo_state_var_evol_rate_name").size()),
  _num_mat_state_var(parameters.get<std::vector<std::string> >("uo_state_var_name").size())
{
  _mat_prop_state_var_local.resize(_num_mat_state_var);
  _mat_prop_state_var_local_old.resize(_num_mat_state_var);
  _mat_prop_state_var_evol_rates.resize(_num_mat_state_var_evol_rates);

  for (unsigned int i = 0 ; i < _num_mat_state_var; ++i)
  {
    _mat_prop_state_var_local[i] = &getMaterialProperty< std::vector<Real> >(parameters.get<std::vector<std::string> >("uo_state_var_name")[i]+"_local");
    _mat_prop_state_var_local_old[i] = &getMaterialProperty< std::vector<Real> >(parameters.get<std::vector<std::string> >("uo_state_var_name")[i]+"_local_old");
  }
  
  for (unsigned int i = 0 ; i < _num_mat_state_var_evol_rates ; ++i)
  {   
    _mat_prop_state_var_evol_rates[i] = &getMaterialProperty< std::vector<Real> >(parameters.get<std::vector<std::string> >("uo_state_var_evol_rate_name")[i]);
  }
}

bool
CrystalPlasticityStateVariableGSS::updateStateVariable(unsigned int qp, Real dt, std::vector<Real> & val) const
{
  
  std::vector<Real> hprops = _uo_slip_resistance->getHardnessParams();

  Real r = hprops[0];
  Real h0 = hprops[1];
  Real tau_init = hprops[2];
  Real tau_sat = hprops[3];

  DenseVector<Real> hb(_nss);
  Real qab;

  val.resize(_nss);

  Real a = hprops[4]; // Kalidindi

  for (unsigned int i = 0; i < _nss; ++i)
    hb(i) = h0 * std::pow(std::abs(1.0 - (*_mat_prop_state_var_local[0])[qp][i]/tau_sat),a) * copysign(1.0,1.0 - (*_mat_prop_state_var_local[0])[qp][i]/tau_sat);

  for (unsigned int i=0; i < _nss; ++i)
  {
    val[i] = (*_mat_prop_state_var_local_old[0])[qp][i];
    for (unsigned int j = 0; j < _nss; ++j)
    {
      unsigned int iplane, jplane;
      iplane = i/3;
      jplane = j/3;

      if (iplane == jplane) // Kalidindi
        qab = 1.0;
      else
        qab = r;

      val[i] += std::abs((*_mat_prop_state_var_evol_rates[0])[qp][j]) * qab * hb(j);
      //std::cout << "val[ " << i << "] ===> " << (*_mat_prop_state_var_local_old[0])[qp][i] << std::endl;
    }
  }
  return true;
}

void
CrystalPlasticityStateVariableGSS::initSlipSysProps(std::vector<Real> & val) const
{
  _uo_slip_resistance->initSlipSysProps(val);
}
