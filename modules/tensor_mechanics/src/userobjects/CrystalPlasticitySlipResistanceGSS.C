/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
//  phenomenological constitutive models' slip resistance userobject class.
//
#include "CrystalPlasticitySlipResistanceGSS.h"

template<>
InputParameters validParams<CrystalPlasticitySlipResistanceGSS>()
{
  InputParameters params = validParams<CrystalPlasticitySlipResistance>();
  params.addParam<FileName>("slip_sys_res_prop_file_name", "", "Name of the file containing the initial values of slip system resistances");
  MooseEnum intvar_read_options("slip_sys_file slip_sys_res_file none","none");
  params.addParam<FileName>("slip_sys_hard_prop_file_name", "", "Name of the file containing the values of hardness evolution parameters");
  params.addParam<MooseEnum>("intvar_read_type", intvar_read_options, "Read from options for initial value of internal variables: Default from .i file");
  params.addParam<std::vector<Real> >("gprops", "Initial values of slip system resistances");
  params.addParam<std::vector<Real> >("hprops", "Hardening properties");
  params.addParam<std::vector<std::string> >("uo_state_var_name", "Name of state variable property: Same as state variable user object specified in input file.");
  params.addParam<std::string>("uo_slip_rate_name", "Name of slip rate property: Same as slip rate user object specified in input file.");
  params.addClassDescription("Phenomenological constitutive models' slip resistance base class.  Override the virtual functions in your class");
  return params;
}

CrystalPlasticitySlipResistanceGSS::CrystalPlasticitySlipResistanceGSS(const InputParameters & parameters) :
  CrystalPlasticitySlipResistance(parameters),
  _slip_sys_res_prop_file_name(getParam<FileName>("slip_sys_res_prop_file_name")),
  _slip_sys_hard_prop_file_name(getParam<FileName>("slip_sys_hard_prop_file_name")),
  _intvar_read_type(getParam<MooseEnum>("intvar_read_type")),
  _gprops(getParam<std::vector<Real> >("gprops")),
  _hprops(getParam<std::vector<Real> >("hprops")),
  _num_mat_state_var(parameters.get<std::vector<std::string> >("uo_state_var_name").size()),
  _mat_prop_slip_rate(getMaterialProperty<std::vector<Real> >(parameters.get<std::string>("uo_slip_rate_name")))
{
  _mat_prop_state_var_local.resize(_num_mat_state_var);

  for (unsigned int i = 0 ; i < _num_mat_state_var; ++i)
  {
    _mat_prop_state_var_local[i] = &getMaterialProperty< std::vector<Real> >(parameters.get<std::vector<std::string> >("uo_state_var_name")[i]+"_local");
  }

  if (_slip_sys_hard_prop_file_name.length()!=0)
    readFileHardnessParams();
  else
    assignHardnessParams();
}

void
CrystalPlasticitySlipResistanceGSS::initSlipSysProps(std::vector<Real> & val) const
{
  switch (_intvar_read_type)
  { 
    case 0:
      assignSlipSysRes(val);
      break;
    case 1: 
      readFileInitSlipSysRes(val);
      break;
    default:
      getInitSlipSysRes(val);
  }
}

void
CrystalPlasticitySlipResistanceGSS::readFileInitSlipSysRes(std::vector<Real> & val) const
{
  MooseUtils::checkFileReadable(_slip_sys_res_prop_file_name);

  std::ifstream file;
  file.open(_slip_sys_res_prop_file_name.c_str());

  for (unsigned int i = 0; i < _nss; ++i)
    if (!(file >> val[i]))
      mooseError("Error FiniteStrainUObasedCP: Premature end of slip_sys_res_prop file");

  file.close();
}

void
CrystalPlasticitySlipResistanceGSS::assignSlipSysRes(std::vector<Real> & val) const
{
}

// Read initial slip system resistances  from .i file
void
CrystalPlasticitySlipResistanceGSS::getInitSlipSysRes(std::vector<Real> & val) const
{
  if (_gprops.size() <= 0)
    mooseError("FiniteStrainUObasedCP: Error in reading slip system resistance properties: Specify input in .i file or in slip_sys_res_prop_file or in slip_sys_file");

  val.resize(_nss, 0.0);

  unsigned int num_data_grp = 3; //Number of data per group e.g. start_slip_sys, end_slip_sys, value

  for (unsigned int i = 0; i < _gprops.size()/num_data_grp; ++i)
  {
    Real vs,ve;
    unsigned int is, ie;

    vs = _gprops[i * num_data_grp];
    ve = _gprops[i * num_data_grp + 1];

    if (vs <= 0 || ve <= 0)
      mooseError( "FiniteStrainUObasedCP: Indices in gss property read must be positive integers: is = " << vs << " ie = " << ve );

    if (vs != floor(vs) || ve != floor(ve))
      mooseError("FiniteStrainUObasedCP: Error in reading slip system resistances: Values specifying start and end number of slip system groups should be integer");

    is = static_cast<unsigned int>(vs);
    ie = static_cast<unsigned int>(ve);

    if (is > ie)
      mooseError("FiniteStrainUObasedCP: Start index is = " << is << " should be greater than end index ie = " << ie << " in slip system resistance property read");

    for (unsigned int j = is; j <= ie; ++j)
      val[j-1] = _gprops[i * num_data_grp + 2];
  } 

  for (unsigned int i = 0; i < _nss; ++i)
      if (val[i] <= 0.0)
      mooseError("FiniteStrainUObasedCP: Value of resistance for slip system " << i + 1 << " non positive");
}

void
CrystalPlasticitySlipResistanceGSS::readFileHardnessParams()
{   
} 

void
CrystalPlasticitySlipResistanceGSS::assignHardnessParams()
{
  if (_hprops.size() <= 0)
    mooseError("FiniteStrainUObasedCP: Error in reading hardness properties: Specify input in .i file or a slip_sys_hard_prop_file_name");

  _r = _hprops[0];
  _h0 = _hprops[1];
  _tau_init = _hprops[2];
  _tau_sat = _hprops[3];
}

std::vector<Real>
CrystalPlasticitySlipResistanceGSS::getHardnessParams() const
{
  return _hprops;
}

bool
CrystalPlasticitySlipResistanceGSS::calcSlipResistance(unsigned int qp, std::vector<Real> & val) const
{
  val.resize(_nss);
  for (unsigned int i = 0; i < _nss; i++)
    val[i] = 0.0;
  return true;
}
