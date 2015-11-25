/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "FiniteStrainUObasedCP.h"
#include "petscblaslapack.h"

#include "CrystalPlasticitySlipRate.h"
#include "CrystalPlasticitySlipResistance.h"
#include "CrystalPlasticityStateVariable.h"
#include "CrystalPlasticityStateVariableEvolutionRate.h"
#include "CrystalPlasticityStateVariableEvolutionRateComponent.h"

template<>
InputParameters validParams<FiniteStrainUObasedCP>()
{
  InputParameters params = validParams<FiniteStrainMaterial>();
  params.addClassDescription("Crystal Plasticity base class: FCC system with power law flow rule implemented");
  params.addRequiredParam<int >("nss", "Number of slip systems");
  params.addParam<Real>("rtol", 1e-6, "Constitutive stress residue relative tolerance");
  params.addParam<Real>("abs_tol", 1e-6, "Constitutive stress residue absolute tolerance");
  params.addParam<Real>("gtol", 1e2, "Constitutive slip system resistance residual tolerance");
  params.addParam<Real>("slip_incr_tol", 2e-2, "Maximum allowable slip in an increment");
  params.addParam<Real>("zero_tol", 1e-12, "Tolerance for residual check when variable value is zero");
  params.addParam<unsigned int>("maxiter", 100 , "Maximum number of iterations for stress update");
  params.addParam<unsigned int>("maxitergss", 100 , "Maximum number of iterations for slip system resistance update");
  params.addParam<UserObjectName>("read_prop_user_object","The ElementReadPropertyFile GeneralUserObject to read element specific property values from file");
  MooseEnum tan_mod_options("exact none","none");//Type of read
  params.addParam<MooseEnum>("tan_mod_type", tan_mod_options, "Type of tangent moduli for preconditioner: default elastic");
  params.addParam<bool>("save_euler_angle", false , "Saves the Euler angles as Material Property if true");
  params.addParam<bool>("gen_random_stress_flag", false, "Flag to generate random stress to perform time cutback on constitutive failure");
  params.addParam<bool>("input_random_scaling_var", false, "Flag to input scaling variable: _Cijkl(0,0,0,0) when false");
  params.addParam<Real>("random_scaling_var", 1e9, "Random scaling variable: Large value can cause non-positive definiteness");
  params.addParam<unsigned int>("random_seed", 2000, "Random integer used to generate random stress when constitutive failure occurs");
  params.addParam<unsigned int>("maximum_substep_iteration", 1, "Maximum number of substep iteration");
  params.addParam<bool>("use_line_search", false, "Use line search in constitutive update");
  params.addParam<Real>("min_line_search_step_size", 0.01, "Minimum line search step size");
  params.addParam<Real>("line_search_tol",0.5,"Line search bisection method tolerance");
  params.addParam<unsigned int>("line_search_maxiter",20,"Line search bisection method maximum number of iteration");
  MooseEnum line_search_method("CUT_HALF BISECTION","CUT_HALF");
  params.addParam<MooseEnum>("line_search_method",line_search_method,"The method used in line search");
  params.addRequiredParam<std::vector<UserObjectName> >("uo_slip_rates", "List of names of user objects that define the slip rates for this material.");
  params.addRequiredParam<std::vector<UserObjectName> >("uo_slip_resistances", "List of names of user objects that define the slip resistances for this material.");
  params.addRequiredParam<std::vector<UserObjectName> >("uo_state_vars", "List of names of user objects that define the state variable for this material.");
  params.addRequiredParam<std::vector<UserObjectName> >("uo_state_var_evol_rates", "List of names of user objects that define the state variable evolution rates for this material.");
  params.addRequiredParam<std::vector<UserObjectName> >("uo_state_var_evol_rate_comps", "List of names of user objects that define the state variable evolution rate components for this material.");
  return params;
}

FiniteStrainUObasedCP::FiniteStrainUObasedCP(const InputParameters & parameters) :
    FiniteStrainMaterial(parameters),
    _num_uo_slip_rates(parameters.get<std::vector<UserObjectName> >("uo_slip_rates").size()),
    _num_uo_slip_resistances(parameters.get<std::vector<UserObjectName> >("uo_slip_resistances").size()),
    _num_uo_state_vars(parameters.get<std::vector<UserObjectName> >("uo_state_vars").size()),
    _num_uo_state_var_evol_rates(parameters.get<std::vector<UserObjectName> >("uo_state_var_evol_rates").size()),
    _num_uo_state_var_evol_rate_comps(parameters.get<std::vector<UserObjectName> >("uo_state_var_evol_rate_comps").size()),
    _nss(getParam<int>("nss")),
    _rtol(getParam<Real>("rtol")),
    _abs_tol(getParam<Real>("abs_tol")),
    _gtol(getParam<Real>("gtol")),
    _slip_incr_tol(getParam<Real>("slip_incr_tol")),
    _zero_tol(getParam<Real>("zero_tol")),
    _maxiter(getParam<unsigned int>("maxiter")),
    _maxiterg(getParam<unsigned int>("maxitergss")),
    _read_prop_user_object(isParamValid("read_prop_user_object") ? & getUserObject<ElementPropertyReadFile>("read_prop_user_object") : NULL),
    _tan_mod_type(getParam<MooseEnum>("tan_mod_type")),
    _save_euler_angle(getParam<bool>("save_euler_angle")),
    _gen_rndm_stress_flag(getParam<bool>("gen_random_stress_flag")),
    _input_rndm_scale_var(getParam<bool>("input_random_scaling_var")),
    _rndm_scale_var(getParam<Real>("random_scaling_var")),
    _rndm_seed(getParam<unsigned int>("random_seed")),
    _max_substep_iter(getParam<unsigned int>("maximum_substep_iteration")),
    _use_line_search(getParam<bool>("use_line_search")),
    _min_lsrch_step(getParam<Real>("min_line_search_step_size")),
    _lsrch_tol(getParam<Real>("line_search_tol")),
    _lsrch_max_iter(getParam<unsigned int>("line_search_maxiter")),
    _lsrch_method(getParam<MooseEnum>("line_search_method")),
    _fp(declareProperty<RankTwoTensor>("fp")), // Plastic deformation gradient
    _fp_old(declarePropertyOld<RankTwoTensor>("fp")), // Plastic deformation gradient of previous increment
    _pk2(declareProperty<RankTwoTensor>("pk2")), // 2nd Piola Kirchoff Stress
    _pk2_old(declarePropertyOld<RankTwoTensor>("pk2")), // 2nd Piola Kirchoff Stress of previous increment
    _pk2_tmp(declareProperty<RankTwoTensor>("pk2_tmp")),
    _pk2_tmp_old(declareProperty<RankTwoTensor>("pk2_tmp_old")),
    _lag_e(declareProperty<RankTwoTensor>("lage")), // Lagrangian strain
    _lag_e_old(declarePropertyOld<RankTwoTensor>("lage")), // Lagrangian strain of previous increment
    _acc_slip(declareProperty<Real>("acc_slip")), // Accumulated slip
    _acc_slip_old(declarePropertyOld<Real>("acc_slip")), // Accumulated alip of previous increment
    _update_rot(declareProperty<RankTwoTensor>("update_rot")), // Rotation tensor considering material rotation and crystal orientation
    _update_rot_old(declarePropertyOld<RankTwoTensor>("update_rot")),
    _deformation_gradient_old(declarePropertyOld<RankTwoTensor>("deformation_gradient")),
    _slip_incr(_nss),
    _dslipdtau(_nss)
{
  if (_save_euler_angle)
  {
    _euler_ang = &declareProperty< std::vector<Real> >("euler_ang");
    _euler_ang_old = &declarePropertyOld< std::vector<Real> >("euler_ang");
  }

  if (!_input_rndm_scale_var)
    _rndm_scale_var = _Cijkl(0,0,0,0);

  _err_tol = false;

  _delta_dfgrd.zero();

  _first_step_iter = false;
  _last_step_iter = false;
  //Initialize variables in the first iteration of substepping
  _first_substep = true;

  RankTwoTensor::initRandom( _rndm_seed );

  //resize the material properties for each userobject
  _mat_prop_slip_rates.resize(_num_uo_slip_rates);
  _mat_prop_slip_resistances.resize(_num_uo_slip_resistances);
  _mat_prop_state_vars.resize(_num_uo_state_vars);
  _mat_prop_state_vars_old.resize(_num_uo_state_vars);
  _mat_prop_state_var_evol_rates.resize(_num_uo_state_var_evol_rates);
  _mat_prop_state_var_evol_rate_comps.resize(_num_uo_state_var_evol_rate_comps);
  
  //resize the schmid tensor
  _s0.resize(_num_uo_slip_rates);

  for (unsigned int i = 0; i < _num_uo_slip_rates; i++)
    _s0[i].resize(_nss);

  //resize local state variables
  _state_var.resize(_num_uo_state_var_evol_rates);
  _state_var_old.resize(_num_uo_state_var_evol_rates);
  _state_var_prev.resize(_num_uo_state_var_evol_rates);

  _uo_slip_rates.resize(_num_uo_slip_rates);
  _uo_slip_resistances.resize(_num_uo_slip_resistances);
  _uo_state_vars.resize(_num_uo_state_vars);
  _uo_state_var_evol_rates.resize(_num_uo_state_var_evol_rates);
  _uo_state_var_evol_rate_comps.resize(_num_uo_state_var_evol_rate_comps);

  // assign the userobjects
  UserObjectInterface uoi(parameters);
  for (unsigned int i = 0 ; i < _num_uo_slip_rates ; ++i)
  {
    _uo_slip_rates[i] = &uoi.getUserObjectByName<CrystalPlasticitySlipRate>(parameters.get<std::vector<UserObjectName> >("uo_slip_rates")[i]);
    _mat_prop_slip_rates[i] = &declareProperty< std::vector<Real> >(parameters.get<std::vector<UserObjectName> >("uo_slip_rates")[i]);
  }

  for (unsigned int i = 0 ; i < _num_uo_slip_resistances ; ++i)
  {
    _uo_slip_resistances[i] = &uoi.getUserObjectByName<CrystalPlasticitySlipResistance>(parameters.get<std::vector<UserObjectName> >("uo_slip_resistances")[i]);
    _mat_prop_slip_resistances[i] = &declareProperty< std::vector<Real> >(parameters.get<std::vector<UserObjectName> >("uo_slip_resistances")[i]);
  }

  for (unsigned int i = 0 ; i < _num_uo_state_vars ; ++i)
  {
    _uo_state_vars[i] = &uoi.getUserObjectByName<CrystalPlasticityStateVariable>(parameters.get<std::vector<UserObjectName> >("uo_state_vars")[i]);
    _mat_prop_state_vars[i] = &declareProperty< std::vector<Real> >(parameters.get<std::vector<UserObjectName> >("uo_state_vars")[i]);
    _mat_prop_state_vars_old[i] = &declarePropertyOld< std::vector<Real> >(parameters.get<std::vector<UserObjectName> >("uo_state_vars")[i]);

    _state_var[i] = &declareProperty< std::vector<Real> >(parameters.get<std::vector<UserObjectName> >("uo_state_vars")[i] + "_local");
    _state_var_old[i] = &declareProperty< std::vector<Real> >(parameters.get<std::vector<UserObjectName> >("uo_state_vars")[i] + "_local_old");
    _state_var_prev[i] = &declareProperty< std::vector<Real> >(parameters.get<std::vector<UserObjectName> >("uo_state_vars")[i] + "_local_prev");
  }

  for (unsigned int i = 0 ; i < _num_uo_state_var_evol_rates ; ++i)
  {
    _uo_state_var_evol_rates[i] = &uoi.getUserObjectByName<CrystalPlasticityStateVariableEvolutionRate>(parameters.get<std::vector<UserObjectName> >("uo_state_var_evol_rates")[i]);
    _mat_prop_state_var_evol_rates[i] = &declareProperty< std::vector<Real> >(parameters.get<std::vector<UserObjectName> >("uo_state_var_evol_rates")[i]);
  }

  for (unsigned int i = 0 ; i < _num_uo_state_var_evol_rate_comps ; ++i)
  {
    _uo_state_var_evol_rate_comps[i] = &uoi.getUserObjectByName<CrystalPlasticityStateVariableEvolutionRateComponent>(parameters.get<std::vector<UserObjectName> >("uo_state_var_evol_rate_comps")[i]);
    _mat_prop_state_var_evol_rate_comps[i] = &declareProperty< std::vector<Real> >(parameters.get<std::vector<UserObjectName> >("uo_state_var_evol_rate_comps")[i]);
  }
}

void FiniteStrainUObasedCP::initQpStatefulProperties()
{
  for (unsigned int i = 0 ; i < _num_uo_slip_rates ; ++i)
  { 
    (*_mat_prop_slip_rates[i])[_qp].resize(_uo_slip_rates[i]->dataSize());
  }

  for (unsigned int i = 0 ; i < _num_uo_slip_resistances ; ++i)
  {
    (*_mat_prop_slip_resistances[i])[_qp].resize(_uo_slip_resistances[i]->dataSize());
  }

  for (unsigned int i = 0 ; i < _num_uo_state_vars ; ++i)
  {
    (*_mat_prop_state_vars[i])[_qp].resize(_uo_state_vars[i]->dataSize());
    (*_mat_prop_state_vars_old[i])[_qp].resize(_uo_state_vars[i]->dataSize());
  }

  for (unsigned int i = 0 ; i < _num_uo_state_var_evol_rates ; ++i)
  {
    (*_mat_prop_state_var_evol_rates[i])[_qp].resize(_uo_state_var_evol_rates[i]->dataSize());
  }

  for (unsigned int i = 0 ; i < _num_uo_state_var_evol_rate_comps ; ++i)
  {
    (*_mat_prop_state_var_evol_rate_comps[i])[_qp].resize(_uo_state_var_evol_rate_comps[i]->dataSize());
  }

  _stress[_qp].zero();

  _fp[_qp].zero();
  _fp[_qp].addIa(1.0);

  _pk2[_qp].zero();
  _pk2_tmp[_qp].zero();
  _pk2_tmp_old[_qp].zero();

  _acc_slip[_qp] = 0.0;
  _lag_e[_qp].zero();

  _update_rot[_qp].zero();
  _update_rot[_qp].addIa(1.0);

  if (_save_euler_angle)
  {
    (*_euler_ang)[_qp].resize(LIBMESH_DIM);
    (*_euler_ang_old)[_qp].resize(LIBMESH_DIM);
  }

  for (unsigned int i = 0 ; i < _num_uo_state_vars ; ++i)
  {
    std::vector<Real> val;
    _uo_state_vars[i]->initSlipSysProps(val); // Initializes slip system related properties
    (*_mat_prop_state_vars[i])[_qp] = val;
    (*_mat_prop_state_vars_old[i])[_qp] = val;
  }
}

void
FiniteStrainUObasedCP::getEulerAngles()
{
  if (_read_prop_user_object)
  {
    _Euler_angles(0) = _read_prop_user_object->getData(_current_elem, 0);
    _Euler_angles(1) = _read_prop_user_object->getData(_current_elem, 1);
    _Euler_angles(2) = _read_prop_user_object->getData(_current_elem, 2);
  }
}

// Calculate crystal rotation tensor from Euler Angles
void
FiniteStrainUObasedCP::getEulerRotations()
{
  Real phi1, phi, phi2;
  Real cp, cp1, cp2, sp, sp1, sp2;
  RankTwoTensor RT;
  Real pi = libMesh::pi;

  phi1 = _Euler_angles(0) * (pi/180.0);
  phi =  _Euler_angles(1) * (pi/180.0);
  phi2 = _Euler_angles(2) * (pi/180.0);

  cp1 = std::cos(phi1);
  cp2 = std::cos(phi2);
  cp = std::cos(phi);

  sp1 = std::sin(phi1);
  sp2 = std::sin(phi2);
  sp = std::sin(phi);

  RT(0,0) = cp1 * cp2 - sp1 * sp2 * cp;
  RT(0,1) = sp1 * cp2 + cp1 * sp2 * cp;
  RT(0,2) = sp2 * sp;
  RT(1,0) = -cp1 * sp2 - sp1 * cp2 * cp;
  RT(1,1) = -sp1 * sp2 + cp1 * cp2 * cp;
  RT(1,2) = cp2 * sp;
  RT(2,0) = sp1 * sp;
  RT(2,1) = -cp1 * sp;
  RT(2,2) = cp;

  _crysrot = RT.transpose();
}

/**
 * Solves stress residual equation using NR.
 * Updates slip system resistances iteratively.
 */
void FiniteStrainUObasedCP::computeQpStress()
{
  unsigned int substep_iter = 1;//Depth of substepping; Limited to maximum substep iteration
  unsigned int num_substep = 1;//Calculated from substep_iter as 2^substep_iter
  Real dt_original = _dt;//Stores original _dt; Reset at the end of solve
  _first_substep = true;//Initialize variables at substep_iter = 1

  if (_max_substep_iter > 1)
  {
    _dfgrd_tmp_old = _deformation_gradient_old[_qp];
    if (_dfgrd_tmp_old.det() == 0)
      _dfgrd_tmp_old.addIa(1.0);

    _delta_dfgrd = _deformation_gradient[_qp] - _dfgrd_tmp_old;
    _err_tol = true;//Indicator to continue substepping
  }

  //Substepping loop
  while (_err_tol && _max_substep_iter > 1)
  {
    _dt = dt_original/num_substep;

    for (unsigned int istep = 0; istep < num_substep; ++istep)
    {
      _first_step_iter = false;
      if (istep == 0)
        _first_step_iter = true;

      _last_step_iter = false;
      if (istep == num_substep - 1)
        _last_step_iter = true;

      _dfgrd_scale_factor = (static_cast<Real>(istep)+1)/num_substep;

      preSolveQp();
      solveQp();

      if (_err_tol)
      {
        substep_iter++;
        num_substep*=2;
        break;
      }
    }

    _first_substep = false;//Prevents reinitialization
    _dt = dt_original;//Resets dt

#ifdef DEBUG
    if (substep_iter > _max_substep_iter)
      mooseWarning("FiniteStrainUObasedCP: Failure with substepping");
#endif

    if (!_err_tol || substep_iter > _max_substep_iter)
      postSolveQp();//Evaluate variables after successful solve or indicate failure
  }

  //No substepping
  if (_max_substep_iter == 1)
  {
    preSolveQp();
    solveQp();
    postSolveQp();
  }
}

void
FiniteStrainUObasedCP::preSolveQp()
{
  //Initialize variable
  if (_first_substep)
  {
    _Jacobian_mult[_qp].zero();//Initializes jacobian for preconditioner
    getEulerAngles();
    getEulerRotations();

    //calc_schmid_tensor();
    for (unsigned int i = 0 ; i < _num_uo_slip_rates ; ++i)
    {
      _uo_slip_rates[i]->calcSchmidTensor(_crysrot, _s0[i]);
    }

    RealTensorValue rot;

    for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
      for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
        rot(i,j) = _crysrot(i,j);

    _elasticity_tensor[_qp] = _Cijkl;
    _elasticity_tensor[_qp].rotate(rot);
  }

  if (_max_substep_iter == 1)
    _dfgrd_tmp = _deformation_gradient[_qp];//Without substepping
  else
    _dfgrd_tmp = _dfgrd_scale_factor * _delta_dfgrd + _dfgrd_tmp_old;

  _err_tol = false;
}

void
FiniteStrainUObasedCP::solveQp()
{
  preSolveStatevar();
  solveStatevar();
  if (_err_tol)
    return;
  postSolveStatevar();
}

void
FiniteStrainUObasedCP::postSolveQp()
{
  if (_err_tol)
  {
    _err_tol = false;
    if ( _gen_rndm_stress_flag )
      _stress[_qp] = RankTwoTensor::genRandomSymmTensor( _rndm_scale_var, 1.0 );
    else
      mooseError("FiniteStrainUObasedCP: Constitutive failure");
  }
  else
  {
    _stress[_qp] = _fe * _pk2[_qp] * _fe.transpose()/_fe.det();

    _Jacobian_mult[_qp] += calcTangentModuli();//Calculate jacobian for preconditioner

    RankTwoTensor iden;
    iden.addIa(1.0);

    _lag_e[_qp] = _deformation_gradient[_qp].transpose() * _deformation_gradient[_qp] - iden;
    _lag_e[_qp] = _lag_e[_qp] * 0.5;

    RankTwoTensor rot;
    rot = get_current_rotation(_deformation_gradient[_qp]); // Calculate material rotation
    _update_rot[_qp] = rot * _crysrot;

    if (_save_euler_angle)
      for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
        (*_euler_ang)[_qp][i] = _Euler_angles(i);
  }
}

void
FiniteStrainUObasedCP::preSolveStatevar()
{
  if (_max_substep_iter == 1)
    _first_step_iter = true;

  if (_first_step_iter)
  {
    for (unsigned int i = 0; i < _num_uo_state_vars; i++)
    {
      (*_state_var[i])[_qp] = (*_state_var_old[i])[_qp] = (*_mat_prop_state_vars_old[i])[_qp];
    }
  }
  else
  {
    for (unsigned int i = 0; i < _num_uo_state_vars; i++)
    {
      (*_state_var[i])[_qp] = (*_state_var_old[i])[_qp];
    }
  }

  for (unsigned int i = 0; i < _num_uo_slip_rates; i++)
  {
    std::vector<Real> val;
    _uo_slip_resistances[i]->calcSlipResistance(_qp, val);
    (*_mat_prop_slip_resistances[i])[_qp] = val;
  }
}

void
FiniteStrainUObasedCP::solveStatevar()
{
  Real gmax, gdiff;
  unsigned int iterg;
  bool iter_flag = true;
  gmax = 1.1 * _gtol;
  
  iterg = 0;

  while (iter_flag && iterg < _maxiterg) // Check for slip system resistance update tolerance
  {
    preSolveStress();
    solveStress();
    if (_err_tol)
      return;
    postSolveStress();

    update_slip_system_resistance(); // Update slip system resistance

    for (unsigned int i = 0; i < _num_uo_state_vars; i++)
    {
      iter_flag = getIterFlagVar((*_state_var[i])[_qp], (*_state_var_prev[i])[_qp], (*_state_var_old[i])[_qp]);
      if (iter_flag)
        break;
    }
    
    iterg++;
  }

  if (iterg == _maxiterg)
  {
#ifdef DEBUG
    mooseWarning("FiniteStrainUObasedCP: Hardness Integration error gmax" << gmax << "\n");
#endif
    _err_tol = true;
  }
}

bool
FiniteStrainUObasedCP::getIterFlagVar(std::vector<Real> var, std::vector<Real> var_prev, std::vector<Real> var_old)
{
  bool iter_flag = false;
  Real gdiff;
  unsigned int n = var.size();

  //Real gmax = 0.0;

  for (unsigned i = 0; i < n; ++i)
  {
    gdiff = std::abs(var[i] - var_prev[i]);//Calculate increment size

    if (std::abs(var_old[i]) < _zero_tol && gdiff > _zero_tol)
      return iter_flag = true;

    if (std::abs(var_old[i]) >  _zero_tol && gdiff > _gtol * std::abs(var_old[i]))
      return iter_flag = true;

    //if ( gdiff > gmax)
    //  gmax = gdiff;
  }

  //if (gmax > _gtol)
  //  return iter_flag = true;

  return iter_flag;
}

void
FiniteStrainUObasedCP::postSolveStatevar()
{
  if (_max_substep_iter == 1)
    _last_step_iter = true;

  if (_last_step_iter)
  {
    for (unsigned int i = 0; i < _num_uo_state_vars; i++)
    {
      (*_mat_prop_state_vars[i])[_qp] = (*_state_var[i])[_qp];
    }
  }
  else
  {
    for (unsigned int i = 0; i < _num_uo_state_vars; i++)
    {
      (*_state_var_old[i])[_qp] = (*_state_var[i])[_qp];
    }
  }
}

void
FiniteStrainUObasedCP::preSolveStress()
{
  if (_max_substep_iter == 1)//No substepping
  {
    _pk2_tmp[_qp] = _pk2_old[_qp];
    _fp_old_inv = _fp_old[_qp].inverse();
    _fp_inv = _fp_old_inv;
    _fp_prev_inv = _fp_inv;
  }
  else
  {
    if (_first_step_iter)
    {
      _pk2_tmp[_qp] = _pk2_tmp_old[_qp] = _pk2_old[_qp];
      _fp_old_inv = _fp_old[_qp].inverse();
    }
    else
      _pk2_tmp[_qp] = _pk2_tmp_old[_qp];

    _fp_inv = _fp_old_inv;
    _fp_prev_inv = _fp_inv;
  }
}

void
FiniteStrainUObasedCP::solveStress()
{
  unsigned int iter = 0;
  RankTwoTensor resid, dpk2;
  RankFourTensor jac;
  Real rnorm, rnorm0, rnorm_prev;

  calc_resid_jacob(resid, jac); // Calculate stress residual
  if (_err_tol)
  {
#ifdef DEBUG
    mooseWarning("FiniteStrainUObasedCP: Slip increment exceeds tolerance - Element number " << _current_elem->id() << " Gauss point = " << _qp);
#endif
    return;
  }

  rnorm = resid.L2norm();
  rnorm0 = rnorm;

  while (rnorm > _rtol * rnorm0 && rnorm0 > _abs_tol && iter <  _maxiter) // Check for stress residual tolerance
  {
    dpk2 = - jac.invSymm() * resid; // Calculate stress increment
    _pk2_tmp[_qp] = _pk2_tmp[_qp] + dpk2; // Update stress
    calc_resid_jacob(resid,jac);
    internalVariableUpdateNRiteration(); //update _fp_prev_inv

    if (_err_tol)
    {
#ifdef DEBUG
      mooseWarning("FiniteStrainUObasedCP: Slip increment exceeds tolerance - Element number " << _current_elem->id() << " Gauss point = " << _qp);
#endif
      return;
    }

    rnorm_prev = rnorm;
    rnorm = resid.L2norm();

    if (_use_line_search && rnorm > rnorm_prev && !line_search_update(rnorm_prev, dpk2))
    {
#ifdef DEBUG
      mooseWarning("FiniteStrainUObasedCP: Failed with line search");
#endif
      _err_tol = true;
      return;
    }

    if (_use_line_search)
      rnorm = resid.L2norm();

    iter++;
  }

  if (iter >= _maxiter)
  {
#ifdef DEBUG
    mooseWarning("FiniteStrainUObasedCP: Stress Integration error rmax = " << rnorm);
#endif
    _err_tol = true;
  }
}

void
FiniteStrainUObasedCP::postSolveStress()
{
  if (_max_substep_iter == 1)//No substepping
  {
    _fp[_qp] = _fp_inv.inverse();
    _pk2[_qp] = _pk2_tmp[_qp];
  }
  else
  {
    if (_last_step_iter)
    {
      _fp[_qp] = _fp_inv.inverse();
      _pk2[_qp] = _pk2_tmp[_qp];
    }
    else
    {
      _fp_old_inv = _fp_inv;
      _pk2_tmp_old[_qp] = _pk2_tmp[_qp];
    }
  }
}

// Update slip system resistance. Overide to incorporate new slip system resistance laws
void
FiniteStrainUObasedCP::update_slip_system_resistance()
{
  for (unsigned int i = 0; i < _num_uo_state_vars; i++)
  {
    (*_state_var_prev[i])[_qp] = (*_state_var[i])[_qp];
  }

  for (unsigned int i = 0; i < _num_uo_state_var_evol_rate_comps; i++)
  {
    std::vector<Real> val;
    _uo_state_var_evol_rate_comps[i]->calcStateVariableEvolutionRateComponent(_qp, val);
    (*_mat_prop_state_var_evol_rate_comps[i])[_qp] = val;
  }

  for (unsigned int i = 0; i < _num_uo_state_var_evol_rates; i++)
  {
    std::vector<Real> val;
    _uo_state_var_evol_rates[i]->calcStateVariableEvolutionRate(_qp, val);
    (*_mat_prop_state_var_evol_rates[i])[_qp] = val;
  }
  
  for (unsigned int i = 0; i < _num_uo_state_vars; i++)
  {
    std::vector<Real> val;
    _uo_state_vars[i]->updateStateVariable(_qp, _dt, val);
    (*_state_var[i])[_qp] = val;
  }

  for (unsigned int i = 0; i < _num_uo_slip_rates; i++)
  {
    std::vector<Real> val;
    _uo_slip_resistances[i]->calcSlipResistance(_qp, val);
    (*_mat_prop_slip_resistances[i])[_qp] = val;
  }
}

// Calculates stress residual equation and jacobian
void
FiniteStrainUObasedCP::calc_resid_jacob( RankTwoTensor & resid, RankFourTensor & jac)
{
  calcResidual( resid );
  if (_err_tol)
    return;
  calcJacobian( jac );
}

void
FiniteStrainUObasedCP::calcResidual( RankTwoTensor &resid )
{
  RankTwoTensor iden, ce, ee, ce_pk2, eqv_slip_incr, pk2_new;

  iden.zero();
  iden.addIa(1.0);

  _fe = _dfgrd_tmp * _fp_prev_inv; // _fp_inv  ==> _fp_prev_inv

  // Calculate dslip
  _slip_incr.zero();
  _dslipdtau.zero();

  eqv_slip_incr.zero();

  for (unsigned int i = 0; i < _num_uo_slip_rates; i++)
  {
    std::vector<Real> val1;
    std::vector<Real> val2;
    _uo_slip_rates[i]->calcSlipRate(_qp, _dt, _s0[i], val1);
    _uo_slip_rates[i]->calcSlipRateDerivative(_qp, _dt, _s0[i], val2);
    for (unsigned int j = 0; j < _nss; j++)
    {
      _slip_incr(j) += val1[j];
      _dslipdtau(j) += val2[j];
      eqv_slip_incr += _s0[i][j] * val1[j];
    }
    (*_mat_prop_slip_rates[i])[_qp] = val1;
  }

  if (_err_tol)
    return;

  eqv_slip_incr = iden - eqv_slip_incr;
  _fp_inv = _fp_old_inv * eqv_slip_incr;
  _fe = _dfgrd_tmp * _fp_inv;

  ce = _fe.transpose() * _fe;
  ee = ce - iden;
  ee *= 0.5;

  pk2_new = _elasticity_tensor[_qp] * ee;

  resid = _pk2_tmp[_qp] - pk2_new;
}

void
FiniteStrainUObasedCP::calcJacobian( RankFourTensor &jac )
{
  RankFourTensor dfedfpinv, deedfe, dfpinvdpk2;

  std::vector<RankTwoTensor> dtaudpk2(_nss), dfpinvdslip(_nss);

  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
        dfedfpinv(i,j,k,j) = _dfgrd_tmp(i,k);

  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
      {
        deedfe(i,j,k,i) = deedfe(i,j,k,i) + _fe(k,j) * 0.5;
        deedfe(i,j,k,j) = deedfe(i,j,k,j) + _fe(k,i) * 0.5;
      }


  for (unsigned int i = 0; i < _num_uo_slip_rates; i++)
  {
    for (unsigned int j = 0; j < _nss; j++)
    {
      dtaudpk2[j] = _s0[i][j];
      dfpinvdslip[j] = - _fp_old_inv * _s0[i][j];
    }

    for (unsigned int j = 0; j < _nss; j++)
      dfpinvdpk2 += (dfpinvdslip[j] * _dslipdtau(j)).outerProduct(dtaudpk2[j]);
  }

  jac = RankFourTensor::IdentityFour() - (_elasticity_tensor[_qp] * deedfe * dfedfpinv * dfpinvdpk2);
}

// Calls getMatRot to perform RU factorization of a tensor.
RankTwoTensor
FiniteStrainUObasedCP::get_current_rotation(const RankTwoTensor & a)
{
  return getMatRot(a);
}

// Performs RU factorization of a tensor
RankTwoTensor
FiniteStrainUObasedCP::getMatRot(const RankTwoTensor & a)
{
  RankTwoTensor rot;
  RankTwoTensor c, diag, evec;
  PetscScalar cmat[LIBMESH_DIM][LIBMESH_DIM], work[10];
  PetscReal w[LIBMESH_DIM];
  PetscBLASInt nd = LIBMESH_DIM,
               lwork = 10,
               info;

  c = a.transpose() * a;

  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      cmat[i][j] = c(i,j);

  LAPACKsyev_("V", "U", &nd, &cmat[0][0], &nd, w, work, &lwork, &info);

  if (info != 0)
    mooseError("FiniteStrainUObasedCP: DSYEV function call in getMatRot function failed");

  diag.zero();

  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    diag(i,i) = std::pow(w[i], 0.5);

  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      evec(i,j) = cmat[i][j];

  rot = a * ((evec.transpose() * diag * evec).inverse());

  return rot;
}

// Calculates tangent moduli which is used for global solve
void
FiniteStrainUObasedCP::computeQpElasticityTensor()
{

}

ElasticityTensorR4
FiniteStrainUObasedCP::calcTangentModuli()
{
  ElasticityTensorR4 tan_mod;

  switch ( _tan_mod_type )
  {
    case 0:
      tan_mod = elastoPlasticTangentModuli();
      break;
    default:
      tan_mod = elasticTangentModuli();
  }

  return tan_mod;
}

ElasticityTensorR4
FiniteStrainUObasedCP::elastoPlasticTangentModuli()
{
  ElasticityTensorR4 tan_mod;
  RankTwoTensor pk2fet, fepk2;
  RankFourTensor deedfe, dsigdpk2dfe;

  // Fill in the matrix stiffness material property

  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
      {
        deedfe(i,j,k,i) = deedfe(i,j,k,i) + _fe(k,j) * 0.5;
        deedfe(i,j,k,j) = deedfe(i,j,k,j) + _fe(k,i) * 0.5;
      }

  dsigdpk2dfe = _fe.mixedProductIkJl(_fe) * _elasticity_tensor[_qp] * deedfe;

  pk2fet = _pk2_tmp[_qp] * _fe.transpose();
  fepk2 = _fe * _pk2_tmp[_qp];

  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      for (unsigned int l = 0; l < LIBMESH_DIM; ++l)
      {
        tan_mod(i,j,i,l) = tan_mod(i,j,i,l) + pk2fet(l,j);
        tan_mod(i,j,j,l) = tan_mod(i,j,j,l) + fepk2(i,l);
      }

  tan_mod += dsigdpk2dfe;

  Real je = _fe.det();
  if (je > 0.0)
    tan_mod /= je;

  return tan_mod;
}


ElasticityTensorR4
FiniteStrainUObasedCP::elasticTangentModuli()
{
  return _elasticity_tensor[_qp];//update jacobian_mult
}

bool
FiniteStrainUObasedCP::line_search_update(const Real rnorm_prev, const RankTwoTensor dpk2)
{
  if (_lsrch_method == "CUT_HALF")
  {
    Real rnorm;
    RankTwoTensor resid;
    Real step = 1.0;

    do
    {
      _pk2_tmp[_qp] = _pk2_tmp[_qp] - step * dpk2;
      step /= 2.0;
      _pk2_tmp[_qp] = _pk2_tmp[_qp] + step * dpk2;

      calcResidual(resid);
      rnorm = resid.L2norm();
    }
    while (rnorm > rnorm_prev && step > _min_lsrch_step);

    if (rnorm > rnorm_prev && step <= _min_lsrch_step)
      return false;

    return true;
  }
  else if (_lsrch_method == "BISECTION")
  {
    unsigned int count = 0;
    Real step_a = 0.0;
    Real step_b = 1.0;
    Real step = 1.0;
    Real s_m = 1000.0;
    Real rnorm = 1000.0;

    RankTwoTensor resid;
    calcResidual(resid);
    Real s_b = resid.doubleContraction(dpk2);
    Real rnorm1 = resid.L2norm();
    _pk2_tmp[_qp] = _pk2_tmp[_qp] - dpk2;
    calcResidual(resid);
    Real s_a = resid.doubleContraction(dpk2);
    Real rnorm0 = resid.L2norm();
    _pk2_tmp[_qp] = _pk2_tmp[_qp] + dpk2;

    if ((rnorm1/rnorm0) < _lsrch_tol || s_a*s_b > 0){
      calcResidual(resid);
      return true;
    }

    while ((rnorm/rnorm0) > _lsrch_tol && count < _lsrch_max_iter)
    {
      _pk2_tmp[_qp] = _pk2_tmp[_qp] - step*dpk2;
      step = 0.5 * (step_b + step_a);
      _pk2_tmp[_qp] = _pk2_tmp[_qp] + step*dpk2;
      calcResidual(resid);
      s_m = resid.doubleContraction(dpk2);
      rnorm = resid.L2norm();

      if (s_m*s_a < 0.0){
        step_b = step;
        s_b = s_m;
      }
      if (s_m*s_b < 0.0){
        step_a = step;
        s_a = s_m;
      }
      count++;
    }

    if ((rnorm/rnorm0) < _lsrch_tol && count < _lsrch_max_iter)
      return true;

    return false;
  }
  else{
    mooseError("Line search meothod is not provided.");
    return false;
  }
}

void
FiniteStrainUObasedCP::internalVariableUpdateNRiteration()
{
  _fp_prev_inv = _fp_inv; // update _fp_prev_inv
}
