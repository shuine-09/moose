/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "FiniteStrainUObasedPolyCP.h"
#include "petscblaslapack.h"

#include "CrystalPlasticitySlipRate.h"
#include "CrystalPlasticitySlipResistance.h"
#include "CrystalPlasticityStateVariable.h"
#include "CrystalPlasticityStateVarRateComponent.h"

#include "GrainTracker.h"

template<>
InputParameters validParams<FiniteStrainUObasedPolyCP>()
{
  InputParameters params = validParams<ComputeStressBase>();
  params.addClassDescription("Crystal Plasticity base class: FCC system with power law flow rule implemented");
  params.addParam<Real>("rtol", 1e-6, "Constitutive stress residue relative tolerance");
  params.addParam<Real>("abs_tol", 1e-6, "Constitutive stress residue absolute tolerance");
  params.addParam<Real>("stol", 1e-2, "Constitutive slip system resistance relative residual tolerance");
  params.addParam<Real>("zero_tol", 1e-12, "Tolerance for residual check when variable value is zero");
  params.addParam<unsigned int>("maxiter", 100 , "Maximum number of iterations for stress update");
  params.addParam<unsigned int>("maxiter_state_variable", 100 , "Maximum number of iterations for state variable update");
  MooseEnum tan_mod_options("exact none","none");// Type of read
  params.addParam<MooseEnum>("tan_mod_type", tan_mod_options, "Type of tangent moduli for preconditioner: default elastic");
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
  params.addRequiredParam<std::vector<UserObjectName> >("uo_state_var_evol_rate_comps", "List of names iof user objects that define the state variable evolution rate components for this material.");
  params.addRequiredCoupledVarWithAutoBuild("v", "var_name_base", "op_num", "Array of coupled variables");
  params.addRequiredParam<UserObjectName>("graintracker_object", "The GrainTracker UserObject to get values from."); 
  return params;
}

FiniteStrainUObasedPolyCP::FiniteStrainUObasedPolyCP(const InputParameters & parameters) :
    ComputeStressBase(parameters),
    _num_uo_slip_rates(parameters.get<std::vector<UserObjectName> >("uo_slip_rates").size()),
    _num_uo_slip_resistances(parameters.get<std::vector<UserObjectName> >("uo_slip_resistances").size()),
    _num_uo_state_vars(parameters.get<std::vector<UserObjectName> >("uo_state_vars").size()),
    _num_uo_state_var_evol_rate_comps(parameters.get<std::vector<UserObjectName> >("uo_state_var_evol_rate_comps").size()),
    _rtol(getParam<Real>("rtol")),
    _abs_tol(getParam<Real>("abs_tol")),
    _stol(getParam<Real>("stol")),
    _zero_tol(getParam<Real>("zero_tol")),
    _maxiter(getParam<unsigned int>("maxiter")),
    _maxiterg(getParam<unsigned int>("maxiter_state_variable")),
    _tan_mod_type(getParam<MooseEnum>("tan_mod_type")),
    _max_substep_iter(getParam<unsigned int>("maximum_substep_iteration")),
    _use_line_search(getParam<bool>("use_line_search")),
    _min_lsrch_step(getParam<Real>("min_line_search_step_size")),
    _lsrch_tol(getParam<Real>("line_search_tol")),
    _lsrch_max_iter(getParam<unsigned int>("line_search_maxiter")),
    _lsrch_method(getParam<MooseEnum>("line_search_method")),
    _fp(declareProperty<std::vector<RankTwoTensor> >("fp")), // Plastic deformation gradient
    _fp_old(declarePropertyOld<std::vector<RankTwoTensor> >("fp")), // Plastic deformation gradient of previous increment
    _pk2(declareProperty<std::vector<RankTwoTensor> >("pk2")), // 2nd Piola Kirchoff Stress
    _pk2_old(declarePropertyOld<std::vector<RankTwoTensor> >("pk2")), // 2nd Piola Kirchoff Stress of previous increment
    _update_rot(declareProperty<std::vector<RankTwoTensor> >("update_rot")), // Rotation tensor considering material rotation and crystal orientation
    _update_rot_old(declarePropertyOld<std::vector<RankTwoTensor> >("update_rot")),
    _lag_e(declareProperty<RankTwoTensor>("lage")), // Lagrangian strain
    _deformation_gradient(getMaterialProperty<RankTwoTensor>("deformation_gradient")),
    _deformation_gradient_old(getMaterialPropertyOld<RankTwoTensor>("deformation_gradient")),
    _crysrot(getMaterialProperty<RankTwoTensor>("crysrot")),
    _grain_tracker(getUserObject<GrainTracker>("graintracker_object")),
    _nop(coupledComponents("v")),
    _vals(_nop)
{
  _err_tol = false;

  _delta_dfgrd.zero();

  // resize user objects
  _uo_slip_rates.resize(_num_uo_slip_rates);
  _uo_slip_resistances.resize(_num_uo_slip_resistances);
  _uo_state_vars.resize(_num_uo_state_vars);
  _uo_state_var_evol_rate_comps.resize(_num_uo_state_var_evol_rate_comps);

  // resize the material properties for each userobject
  _mat_prop_slip_rates.resize(_num_uo_slip_rates);
  _mat_prop_slip_resistances.resize(_num_uo_slip_resistances);
  _mat_prop_state_vars.resize(_num_uo_state_vars);
  _mat_prop_state_vars_old.resize(_num_uo_state_vars);
  _mat_prop_state_var_evol_rate_comps.resize(_num_uo_state_var_evol_rate_comps);

  // resize local state variables
  _state_vars_old.resize(_num_uo_state_vars);
  _state_vars_prev.resize(_num_uo_state_vars);

  // resize the flow direction
  _flow_direction.resize(_num_uo_slip_rates);

  // assign the user objects
  for (unsigned int i = 0; i < _num_uo_slip_rates; ++i)
  {
    _uo_slip_rates[i] = &getUserObjectByName<CrystalPlasticitySlipRate>(parameters.get<std::vector<UserObjectName> >("uo_slip_rates")[i]);
    _mat_prop_slip_rates[i] = &declareProperty< std::vector<Real> >(parameters.get<std::vector<UserObjectName> >("uo_slip_rates")[i]);
    _flow_direction[i] = &declareProperty< std::vector<RankTwoTensor> >(parameters.get<std::vector<UserObjectName> >("uo_slip_rates")[i] + "_flow_direction");
  }

  for (unsigned int i = 0; i < _num_uo_slip_resistances; ++i)
  {
    _uo_slip_resistances[i] = &getUserObjectByName<CrystalPlasticitySlipResistance>(parameters.get<std::vector<UserObjectName> >("uo_slip_resistances")[i]);
    _mat_prop_slip_resistances[i] = &declareProperty< std::vector<Real> >(parameters.get<std::vector<UserObjectName> >("uo_slip_resistances")[i]);
  }

  for (unsigned int i = 0; i < _num_uo_state_vars; ++i)
  {
    _uo_state_vars[i] = &getUserObjectByName<CrystalPlasticityStateVariable>(parameters.get<std::vector<UserObjectName> >("uo_state_vars")[i]);
    _mat_prop_state_vars[i] = &declareProperty< std::vector<Real> >(parameters.get<std::vector<UserObjectName> >("uo_state_vars")[i]);
    _mat_prop_state_vars_old[i] = &declarePropertyOld< std::vector<Real> >(parameters.get<std::vector<UserObjectName> >("uo_state_vars")[i]);
  }

  for (unsigned int i = 0; i < _num_uo_state_var_evol_rate_comps; ++i)
  {
    _uo_state_var_evol_rate_comps[i] = &getUserObjectByName<CrystalPlasticityStateVarRateComponent>(parameters.get<std::vector<UserObjectName> >("uo_state_var_evol_rate_comps")[i]);
    _mat_prop_state_var_evol_rate_comps[i] = &declareProperty< std::vector<Real> >(parameters.get<std::vector<UserObjectName> >("uo_state_var_evol_rate_comps")[i]);
  }

  // Initialize variables
  for (unsigned int op = 0; op < _nop; ++op)
    _vals[op] = &coupledValue("v", op);
}

void FiniteStrainUObasedPolyCP::initQpStatefulProperties()
{  
  _active_ops = &(_grain_tracker.getElementalValues(_current_elem->id()));
  
  /// Active number of order parameters
  _n_active_ops = (*_active_ops).size();
  
  if (_n_active_ops)
  {
    _fp[_qp].resize(_n_active_ops);
    _fp_old[_qp].resize(_n_active_ops);
    _pk2[_qp].resize(_n_active_ops);
    _pk2_old[_qp].resize(_n_active_ops);
    _update_rot[_qp].resize(_n_active_ops);
    _update_rot_old[_qp].resize(_n_active_ops);


    for (unsigned int i = 0; i < _num_uo_state_vars; ++i)
    {
      (*_mat_prop_state_vars[i])[_qp].resize(_uo_state_vars[i]->variableSize() * _n_active_ops);
      (*_mat_prop_state_vars_old[i])[_qp].resize(_uo_state_vars[i]->variableSize() * _n_active_ops);
      _state_vars_old[i].resize(_uo_state_vars[i]->variableSize());
      _state_vars_prev[i].resize(_uo_state_vars[i]->variableSize());
    }

    for( unsigned int op = 0; op < _n_active_ops; ++op)
    {
      _fp[_qp][op].zero();
      _fp[_qp][op].addIa(1.0);

      _pk2[_qp][op].zero();
      _update_rot[_qp][op].zero();
      _update_rot[_qp][op].addIa(1.0);
    }

    for (unsigned int i = 0; i < _num_uo_state_vars; ++i)
    {
      for( unsigned int op = 0; op < _n_active_ops; ++op)
        // Initializes slip system related properties
        _uo_state_vars[i]->initSlipSysProps((*_mat_prop_state_vars[i])[_qp], op);
      
      (*_mat_prop_state_vars_old[i])[_qp] = (*_mat_prop_state_vars[i])[_qp];
    }
  }
  _stress[_qp].zero();
  _lag_e[_qp].zero();
}

/**
 * Solves stress residual equation using NR.
 * Updates slip system resistances iteratively.
 */
void FiniteStrainUObasedPolyCP::computeQpStress()
{
  // Userobject based crystal plasticity does not support face/boundary material property calculation.
  if (isBoundaryMaterial())
    return;

  _active_ops = &(_grain_tracker.getElementalValues(_current_elem->id()));
  
  /// Active number of order parameters
  _n_active_ops = (*_active_ops).size();

  _fe.resize(_n_active_ops);
  _fp_old_inv.resize(_n_active_ops);
  _fp_inv.resize(_n_active_ops);

  // resize non-stateful material property
  for (unsigned int i = 0; i < _num_uo_slip_rates; ++i)
  {
    (*_mat_prop_slip_rates[i])[_qp].resize(_uo_slip_rates[i]->variableSize() * _n_active_ops);
    (*_flow_direction[i])[_qp].resize(_uo_slip_rates[i]->variableSize() * _n_active_ops);
  }

  for (unsigned int i = 0; i < _num_uo_slip_resistances; ++i)
    (*_mat_prop_slip_resistances[i])[_qp].resize(_uo_slip_resistances[i]->variableSize() * _n_active_ops);

  for (unsigned int i = 0; i < _num_uo_state_var_evol_rate_comps; ++i)
    (*_mat_prop_state_var_evol_rate_comps[i])[_qp].resize(_uo_state_var_evol_rate_comps[i]->variableSize() * _n_active_ops);
 
  // Depth of substepping; Limited to maximum substep iteration
  unsigned int substep_iter = 1;
  // Calculated from substep_iter as 2^substep_iter
  unsigned int num_substep = 1;
  // Store original _dt; Reset at the end of solve
  Real dt_original = _dt;

  _dfgrd_tmp_old = _deformation_gradient_old[_qp];
  if (_dfgrd_tmp_old.det() == 0)
    _dfgrd_tmp_old.addIa(1.0);

  _delta_dfgrd = _deformation_gradient[_qp] - _dfgrd_tmp_old;

  // Saves the old stateful properties that is modified during sub stepping
  for (unsigned int op = 0; op < _n_active_ops; ++op)
  {
    _grn_index = op;

    for (unsigned int i = 0; i < _num_uo_state_vars; ++i)
      for (unsigned int j = 0; j < _uo_state_vars[i]->variableSize(); ++j)
        _state_vars_old[i][j] = (*_mat_prop_state_vars_old[i])[_qp][_uo_state_vars[i]->variableSize() * _grn_index + j];

    for (unsigned int i = 0; i < _num_uo_slip_rates; ++i)
      _uo_slip_rates[i]->calcFlowDirection(_qp, (*_flow_direction[i])[_qp], _grn_index);
   
    do
    {
      _err_tol = false;

      preSolveQp();

      _dt = dt_original/num_substep;

      for (unsigned int istep = 0; istep < num_substep; ++istep)
      {
        _dfgrd_tmp =  (static_cast<Real>(istep) + 1) / num_substep * _delta_dfgrd + _dfgrd_tmp_old;

        solveQp();

        if (_err_tol)
        {
          substep_iter++;
          num_substep *= 2;
          break;
        }
      }
      if (substep_iter > _max_substep_iter && _err_tol)
        mooseError("FiniteStrainUObasedPolyCP: Constitutive failure");
    }
    while (_err_tol);

    _dt = dt_original;

    postSolveQp();
  }
  
  finalSolveSressQp();
}

void
FiniteStrainUObasedPolyCP::preSolveQp()
{
  for (unsigned int i = 0; i < _num_uo_state_vars; ++i)
    for (unsigned int j = 0; j < _uo_state_vars[i]->variableSize(); ++j)
      (*_mat_prop_state_vars[i])[_qp][_uo_state_vars[i]->variableSize() * _grn_index + j] = (*_mat_prop_state_vars_old[i])[_qp][_uo_state_vars[i]->variableSize() * _grn_index + j] = _state_vars_old[i][j];

  _pk2[_qp][_grn_index] = _pk2_old[_qp][_grn_index];
  _fp_old_inv[_grn_index] = _fp_old[_qp][_grn_index].inverse();
}

void
FiniteStrainUObasedPolyCP::solveQp()
{
  preSolveStatevar();
  solveStatevar();
  if (_err_tol)
    return;
  postSolveStatevar();
}

void
FiniteStrainUObasedPolyCP::finalSolveSressQp()
{
  Real sum_h = 0.0; 

  _stress[_qp].zero();
  
  // Calculate elasticity tensor
  for (unsigned int op = 0; op < _n_active_ops; ++op)
  {
    // Second position contains the order parameter index
    unsigned int op_index = (*_active_ops)[op].second;

    // Interpolation factor for elasticity tensors
    Real h = (1.0 + std::sin(libMesh::pi * ((*_vals[op_index])[_qp] - 0.5)))/2.0;

    // Sum all rotated elasticity tensors
    sum_h += h;
    
    _stress[_qp] += h * _fe[op] * _pk2[_qp][op] * _fe[op].transpose()/_fe[op].det();
  }

  _stress[_qp] /= sum_h;

  // Calculate jacobian for preconditioner
  calcTangentModuli();

  RankTwoTensor iden;
  iden.addIa(1.0);

  _lag_e[_qp] = _deformation_gradient[_qp].transpose() * _deformation_gradient[_qp] - iden;
  _lag_e[_qp] = _lag_e[_qp] * 0.5;

  RankTwoTensor rot;
  // Calculate material rotation
  _deformation_gradient[_qp].getRUDecompositionRotation(rot);
  //_update_rot[_qp] = rot * _crysrot[_qp];
}

void
FiniteStrainUObasedPolyCP::postSolveQp()
{
  // Restores the the old stateful properties after a successful solve
  for (unsigned int i = 0; i < _num_uo_state_vars; ++i)
    for (unsigned int j = 0; j < _uo_state_vars[i]->variableSize(); ++j)
    (*_mat_prop_state_vars_old[i])[_qp][_uo_state_vars[i]->variableSize() * _grn_index + j] = _state_vars_old[i][j];
}

void
FiniteStrainUObasedPolyCP::preSolveStatevar()
{
  for (unsigned int i = 0; i < _num_uo_state_vars; ++i)
    for (unsigned int j = 0; j < _uo_state_vars[i]->variableSize(); ++j)
      (*_mat_prop_state_vars[i])[_qp][_uo_state_vars[i]->variableSize() * _grn_index + j] = (*_mat_prop_state_vars_old[i])[_qp][_uo_state_vars[i]->variableSize() * _grn_index + j];

  for (unsigned int i = 0; i < _num_uo_slip_resistances; ++i)
    for (unsigned int j = 0; j < _uo_slip_resistances[i]->variableSize(); ++j)
      _uo_slip_resistances[i]->calcSlipResistance(_qp, (*_mat_prop_slip_resistances[i])[_qp], _grn_index);

  _fp_inv[_grn_index] = _fp_old_inv[_grn_index];
}

void
FiniteStrainUObasedPolyCP::solveStatevar()
{
  unsigned int iterg;
  bool iter_flag = true;

  iterg = 0;
  // Check for slip system resistance update tolerance
  while (iter_flag && iterg < _maxiterg)
  {
    preSolveStress();
    solveStress();
    if (_err_tol)
      return;
    postSolveStress();

    // Update slip system resistance and state variable
    updateSlipSystemResistanceAndStateVariable();

    if (_err_tol)
      return;

    iter_flag = isStateVariablesConverged();
    iterg++;
  }

  if (iterg == _maxiterg)
  {
#ifdef DEBUG
    mooseWarning("FiniteStrainUObasedPolyCP: Hardness Integration error\n");
#endif
    _err_tol = true;
  }
}

bool
FiniteStrainUObasedPolyCP::isStateVariablesConverged()
{
  Real diff;

  for (unsigned int i = 0; i < _num_uo_state_vars; ++i)
  {
    for (unsigned j = 0; j < _uo_state_vars[i]->variableSize(); ++ j)
    {
      diff = std::abs((*_mat_prop_state_vars[i])[_qp][_uo_state_vars[i]->variableSize() * _grn_index + j] - _state_vars_prev[i][j]);// Calculate increment size
      if (std::abs((*_mat_prop_state_vars_old[i])[_qp][_uo_state_vars[i]->variableSize() * _grn_index + j]) < _zero_tol && diff > _zero_tol)
        return true;
      if (std::abs((*_mat_prop_state_vars_old[i])[_qp][_uo_state_vars[i]->variableSize() * _grn_index + j]) >  _zero_tol && diff > _stol * std::abs((*_mat_prop_state_vars_old[i])[_qp][_uo_state_vars[i]->variableSize() * _grn_index + j]))
        return true;
    }
  }
  return false;
}

void
FiniteStrainUObasedPolyCP::postSolveStatevar()
{
  for (unsigned int i = 0; i < _num_uo_state_vars; ++i)
    for (unsigned int j = 0; j < _uo_state_vars[i]->variableSize(); ++j)
    (*_mat_prop_state_vars_old[i])[_qp][_uo_state_vars[i]->variableSize() * _grn_index + j] = (*_mat_prop_state_vars[i])[_qp][_uo_state_vars[i]->variableSize() * _grn_index + j];

  _fp_old_inv[_grn_index] = _fp_inv[_grn_index];
}

void
FiniteStrainUObasedPolyCP::preSolveStress()
{
}

void
FiniteStrainUObasedPolyCP::solveStress()
{
  unsigned int iter = 0;
  RankTwoTensor dpk2;
  Real rnorm, rnorm0, rnorm_prev;

  // Calculate stress residual
  calcResidJacob();
  if (_err_tol)
  {
#ifdef DEBUG
    mooseWarning("FiniteStrainUObasedPolyCP: Slip increment exceeds tolerance - Element number " << _current_elem->id() << " Gauss point = " << _qp);
#endif
    return;
  }

  rnorm = _resid.L2norm();
  rnorm0 = rnorm;

  // Check for stress residual tolerance
  while (rnorm > _rtol * rnorm0 && rnorm0 > _abs_tol && iter <  _maxiter)
  {
    // Calculate stress increment
    dpk2 = - _jac.invSymm() * _resid;
    _pk2[_qp][_grn_index] = _pk2[_qp][_grn_index] + dpk2;
    calcResidJacob();

    if (_err_tol)
    {
#ifdef DEBUG
      mooseWarning("FiniteStrainUObasedPolyCP: Slip increment exceeds tolerance - Element number " << _current_elem->id() << " Gauss point = " << _qp);
#endif
      return;
    }

    rnorm_prev = rnorm;
    rnorm = _resid.L2norm();

    if (_use_line_search && rnorm > rnorm_prev && !lineSearchUpdate(rnorm_prev, dpk2))
    {
#ifdef DEBUG
      mooseWarning("FiniteStrainUObasedPolyCP: Failed with line search");
#endif
      _err_tol = true;
      return;
    }

    if (_use_line_search)
      rnorm = _resid.L2norm();

    iter++;
  }

  if (iter >= _maxiter)
  {
#ifdef DEBUG
    mooseWarning("FiniteStrainUObasedPolyCP: Stress Integration error rmax = " << rnorm);
#endif
    _err_tol = true;
  }
}

void
FiniteStrainUObasedPolyCP::postSolveStress()
{
  _fp[_qp][_grn_index] = _fp_inv[_grn_index].inverse();
}

void
FiniteStrainUObasedPolyCP::updateSlipSystemResistanceAndStateVariable()
{
  for (unsigned int i = 0; i < _num_uo_state_vars; ++i)
    for (unsigned int j = 0; j < _uo_state_vars[i]->variableSize(); ++j)
      _state_vars_prev[i][j] = (*_mat_prop_state_vars[i])[_qp][_uo_state_vars[i]->variableSize() * _grn_index + j];

  for (unsigned int i = 0; i < _num_uo_state_var_evol_rate_comps; ++i)
      _uo_state_var_evol_rate_comps[i]->calcStateVariableEvolutionRateComponent(_qp, (*_mat_prop_state_var_evol_rate_comps[i])[_qp], _grn_index);

  for (unsigned int i = 0; i < _num_uo_state_vars; ++i)
  {
    if (!_uo_state_vars[i]->updateStateVariable(_qp, _dt, (*_mat_prop_state_vars[i])[_qp], _grn_index))
      _err_tol = true;
  }

  for (unsigned int i = 0; i < _num_uo_slip_resistances; ++i)
    _uo_slip_resistances[i]->calcSlipResistance(_qp, (*_mat_prop_slip_resistances[i])[_qp], _grn_index);
}

// Calculates stress residual equation and jacobian
void
FiniteStrainUObasedPolyCP::calcResidJacob()
{
  calcResidual();
  if (_err_tol)
    return;
  calcJacobian();
}

void
FiniteStrainUObasedPolyCP::getSlipRates()
{
  for (unsigned int i = 0; i < _num_uo_slip_rates; ++i)
  {
    if (!_uo_slip_rates[i]->calcSlipRate(_qp, _dt, (*_mat_prop_slip_rates[i])[_qp], _grn_index))
    {
      _err_tol = true;
      return;
    }
  }
}

void
FiniteStrainUObasedPolyCP::calcResidual()
{
  RankTwoTensor iden, ce, ee, ce_pk2, eqv_slip_incr, pk2_new;

  iden.zero();
  iden.addIa(1.0);

  eqv_slip_incr.zero();

  getSlipRates();

  if (_err_tol)
    return;

  for (unsigned int i = 0; i < _num_uo_slip_rates; ++i)
    for (unsigned int j = 0; j < _uo_slip_rates[i]->variableSize(); ++j)
      eqv_slip_incr += (*_flow_direction[i])[_qp][_uo_slip_rates[i]->variableSize() * _grn_index +j] * (*_mat_prop_slip_rates[i])[_qp][_uo_slip_rates[i]->variableSize() * _grn_index + j] * _dt;

  eqv_slip_incr = iden - eqv_slip_incr;
  _fp_inv[_grn_index] = _fp_old_inv[_grn_index] * eqv_slip_incr;
  _fe[_grn_index] = _dfgrd_tmp * _fp_inv[_grn_index];

  ce = _fe[_grn_index].transpose() * _fe[_grn_index];
  ee = ce - iden;
  ee *= 0.5;

  pk2_new = _elasticity_tensor[_qp] * ee;

  _resid = _pk2[_qp][_grn_index] - pk2_new;
}

void
FiniteStrainUObasedPolyCP::calcJacobian()
{
  RankFourTensor dfedfpinv, deedfe, dfpinvdpk2;

  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
        dfedfpinv(i,j,k,j) = _dfgrd_tmp(i,k);

  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
      {
        deedfe(i,j,k,i) = deedfe(i,j,k,i) + _fe[_grn_index](k,j) * 0.5;
        deedfe(i,j,k,j) = deedfe(i,j,k,j) + _fe[_grn_index](k,i) * 0.5;
      }

  for (unsigned int i = 0; i < _num_uo_slip_rates; ++i)
  {
    unsigned int nss = _uo_slip_rates[i]->variableSize();
    std::vector<RankTwoTensor> dtaudpk2(nss), dfpinvdslip(nss);
    std::vector<Real> dslipdtau;
    dslipdtau.resize(nss);
    _uo_slip_rates[i]->calcSlipRateDerivative(_qp, _dt, dslipdtau, _grn_index);
    for (unsigned int j = 0; j < nss; j++)
    {
      dtaudpk2[j] = (*_flow_direction[i])[_qp][nss * _grn_index + j];
      dfpinvdslip[j] = - _fp_old_inv[_grn_index] * (*_flow_direction[i])[_qp][nss * _grn_index + j];
    }

    for (unsigned int j = 0; j < nss; j++)
      dfpinvdpk2 += (dfpinvdslip[j] * dslipdtau[j] * _dt).outerProduct(dtaudpk2[j]);
  }
  _jac = RankFourTensor::IdentityFour() - (_elasticity_tensor[_qp] * deedfe * dfedfpinv * dfpinvdpk2);
}

void
FiniteStrainUObasedPolyCP::calcTangentModuli()
{
  switch (_tan_mod_type)
  {
    case 0:
      elastoPlasticTangentModuli();
      break;
    default:
      elasticTangentModuli();
  }
}

void
FiniteStrainUObasedPolyCP::elastoPlasticTangentModuli()
{
  _Jacobian_mult[_qp].zero();
  
  Real sum_h = 0.0;

  for (unsigned int op = 0; op < _n_active_ops; ++op)
  {
    // Second position contains the order parameter index
    unsigned int op_index = (*_active_ops)[op].second;

    // Interpolation factor for elasticity tensors
    Real h = (1.0 + std::sin(libMesh::pi * ((*_vals[op_index])[_qp] - 0.5)))/2.0;

    // Sum all rotated elasticity tensors
    sum_h += h;

    RankFourTensor tan_mod;
    RankTwoTensor pk2fet, fepk2;
    RankFourTensor deedfe, dsigdpk2dfe, dfedf;

    // Fill in the matrix stiffness material property
    for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
      for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
        for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
        {
          deedfe(i,j,k,i) = deedfe(i,j,k,i) + _fe[op](k,j) * 0.5;
          deedfe(i,j,k,j) = deedfe(i,j,k,j) + _fe[op](k,i) * 0.5;
        }

    dsigdpk2dfe = _fe[op].mixedProductIkJl(_fe[op]) * _elasticity_tensor[_qp] * deedfe;

    pk2fet = _pk2[_qp][op] * _fe[op].transpose();
    fepk2 = _fe[op] * _pk2[_qp][op];

    for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
      for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
        for (unsigned int l = 0; l < LIBMESH_DIM; ++l)
        {
          tan_mod(i,j,i,l) += pk2fet(l,j);
          tan_mod(i,j,j,l) += fepk2(i,l);
        }

    tan_mod += dsigdpk2dfe;

    Real je = _fe[op].det();
    if (je > 0.0)
      tan_mod /= je;

    for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
      for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
        for (unsigned int l = 0; l < LIBMESH_DIM; ++l)
          dfedf(i,j,i,l) =  _fp_inv[op](l,j);


    _Jacobian_mult[_qp] += tan_mod * dfedf * h;
  }
  _Jacobian_mult[_qp] /= sum_h;
}

void
FiniteStrainUObasedPolyCP::elasticTangentModuli()
{
  // update jacobian_mult
 _Jacobian_mult[_qp] = _elasticity_tensor[_qp];
}

bool
FiniteStrainUObasedPolyCP::lineSearchUpdate(const Real rnorm_prev, const RankTwoTensor dpk2)
{
  switch (_lsrch_method)
  {
    case 0: // CUT_HALF
    {
      Real rnorm;
      Real step = 1.0;

      do
      {
        _pk2[_qp][_grn_index] = _pk2[_qp][_grn_index] - step * dpk2;
        step /= 2.0;
        _pk2[_qp][_grn_index] = _pk2[_qp][_grn_index] + step * dpk2;

        calcResidual();
        rnorm = _resid.L2norm();
      }
      while (rnorm > rnorm_prev && step > _min_lsrch_step);

      // has norm improved or is the step still above minumum search step size?
      return (rnorm <= rnorm_prev || step > _min_lsrch_step);
    }

    case 1: // BISECTION
    {
      unsigned int count = 0;
      Real step_a = 0.0;
      Real step_b = 1.0;
      Real step = 1.0;
      Real s_m = 1000.0;
      Real rnorm = 1000.0;

      calcResidual();
      Real s_b = _resid.doubleContraction(dpk2);
      Real rnorm1 = _resid.L2norm();
      _pk2[_qp][_grn_index] = _pk2[_qp][_grn_index] - dpk2;
      calcResidual();
      Real s_a = _resid.doubleContraction(dpk2);
      Real rnorm0 = _resid.L2norm();
      _pk2[_qp][_grn_index] = _pk2[_qp][_grn_index] + dpk2;

      if ((rnorm1/rnorm0) < _lsrch_tol || s_a*s_b > 0){
        calcResidual();
        return true;
      }

      while ((rnorm/rnorm0) > _lsrch_tol && count < _lsrch_max_iter)
      {
        _pk2[_qp][_grn_index] = _pk2[_qp][_grn_index] - step*dpk2;
        step = 0.5 * (step_b + step_a);
        _pk2[_qp][_grn_index] = _pk2[_qp][_grn_index] + step*dpk2;
        calcResidual();
        s_m = _resid.doubleContraction(dpk2);
        rnorm = _resid.L2norm();

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

      // below tolerance and max iterations?
      return  ((rnorm/rnorm0) < _lsrch_tol && count < _lsrch_max_iter);
    }

    default:
      mooseError("Line search method is not provided.");
  }
}
