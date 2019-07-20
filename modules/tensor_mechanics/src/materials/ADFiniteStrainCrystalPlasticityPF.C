//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADFiniteStrainCrystalPlasticityPF.h"
#include "petscblaslapack.h"
#include "libmesh/utility.h"
#include "MooseException.h"

#include <fstream>
#include <cmath>

registerADMooseObject("TensorMechanicsApp", ADFiniteStrainCrystalPlasticityPF);

defineADValidParams(
    ADFiniteStrainCrystalPlasticityPF,
    ADComputeStressBase,
    params.addClassDescription(
        "Crystal Plasticity base class: FCC system with power law flow rule implemented");
    params.addRequiredParam<int>("nss", "Number of slip systems");
    params.addParam<std::vector<Real>>("gprops", "Initial values of slip system resistances");
    params.addParam<std::vector<Real>>("hprops", "Hardening properties");
    params.addParam<std::vector<Real>>("flowprops", "Parameters used in slip rate equations");
    params.addRequiredParam<FileName>("slip_sys_file_name",
                                      "Name of the file containing the slip system");
    params.addParam<FileName>(
        "slip_sys_res_prop_file_name",
        "",
        "Name of the file containing the initial values of slip system resistances");
    params.addParam<FileName>(
        "slip_sys_flow_prop_file_name",
        "",
        "Name of the file containing the values of slip rate equation parameters");
    params.addParam<FileName>(
        "slip_sys_hard_prop_file_name",
        "",
        "Name of the file containing the values of hardness evolution parameters");
    params.addParam<Real>("rtol", 1e-6, "Constitutive stress residue relative tolerance");
    params.addParam<Real>("abs_tol", 1e-6, "Constitutive stress residue absolute tolerance");
    params.addParam<Real>("gtol", 1e2, "Constitutive slip system resistance residual tolerance");
    params.addParam<Real>("slip_incr_tol", 2e-2, "Maximum allowable slip in an increment");
    params.addParam<unsigned int>("maxiter", 100, "Maximum number of iterations for stress update");
    params.addParam<unsigned int>("maxitergss",
                                  100,
                                  "Maximum number of iterations for slip system resistance update");
    params.addParam<unsigned int>(
        "num_slip_sys_flowrate_props",
        2,
        "Number of flow rate properties for a slip system"); // Used for reading flow rate
                                                             // parameters
    params.addParam<UserObjectName>("read_prop_user_object",
                                    "The ElementReadPropertyFile "
                                    "GeneralUserObject to read element "
                                    "specific property values from file");
    MooseEnum tan_mod_options("exact none", "none"); // Type of read
    params.addParam<MooseEnum>("tan_mod_type",
                               tan_mod_options,
                               "Type of tangent moduli for preconditioner: default elastic");
    MooseEnum intvar_read_options("slip_sys_file slip_sys_res_file none", "none");
    params.addParam<MooseEnum>(
        "intvar_read_type",
        intvar_read_options,
        "Read from options for initial value of internal variables: Default from .i file");
    params.addParam<unsigned int>("num_slip_sys_props",
                                  0,
                                  "Number of slip system specific properties provided in the file "
                                  "containing slip system normals and directions");
    params.addParam<bool>(
        "gen_random_stress_flag",
        false,
        "Flag to generate random stress to perform time cutback on constitutive failure");
    params.addParam<bool>("input_random_scaling_var",
                          false,
                          "Flag to input scaling variable: _Cijkl(0,0,0,0) when false");
    params.addParam<Real>(
        "random_scaling_var",
        1e9,
        "Random scaling variable: Large value can cause non-positive definiteness");
    params.addParam<unsigned int>(
        "random_seed",
        2000,
        "Random integer used to generate random stress when constitutive failure occurs");
    params.addParam<unsigned int>("maximum_substep_iteration",
                                  1,
                                  "Maximum number of substep iteration");
    params.addParam<bool>("use_line_search", false, "Use line search in constitutive update");
    params.addParam<Real>("min_line_search_step_size", 0.01, "Minimum line search step size");
    params.addParam<Real>("line_search_tol", 0.5, "Line search bisection method tolerance");
    params.addParam<unsigned int>("line_search_maxiter",
                                  20,
                                  "Line search bisection method maximum number of iteration");

    params.addRequiredCoupledVar("c", "Order parameter for damage");
    params.addParam<Real>("kdamage", 0.0, "Stiffness of damaged matrix");
    params.addParam<bool>("use_current_history_variable",
                          false,
                          "Use the current value of the history variable.");
    params.addParam<bool>("beta_p", false, "Include effective plastic work driving energy.");
    params.addParam<bool>("beta_e", false, "Include elastic work driving energy.");
    params.addParam<Real>("W0", 0, "plastic work threshold.");
    params.addParam<Real>("Hall_Petch_const", 0, "Hall_Petch_const");
    params.addParam<Real>("grain_size", 1, "Grain size");
    MooseEnum line_search_method("CUT_HALF BISECTION", "CUT_HALF");
    params.addParam<MooseEnum>("line_search_method",
                               line_search_method,
                               "The method used in line search"););

template <ComputeStage compute_stage>
ADFiniteStrainCrystalPlasticityPF<compute_stage>::ADFiniteStrainCrystalPlasticityPF(
    const InputParameters & parameters)
  : ADComputeStressBase<compute_stage>(parameters),
    _nss(getParam<int>("nss")),
    _gprops(getParam<std::vector<Real>>("gprops")),
    _hprops(getParam<std::vector<Real>>("hprops")),
    _flowprops(getParam<std::vector<Real>>("flowprops")),
    _slip_sys_file_name(getParam<FileName>("slip_sys_file_name")),
    _slip_sys_res_prop_file_name(getParam<FileName>("slip_sys_res_prop_file_name")),
    _slip_sys_flow_prop_file_name(getParam<FileName>("slip_sys_flow_prop_file_name")),
    _slip_sys_hard_prop_file_name(getParam<FileName>("slip_sys_hard_prop_file_name")),
    _rtol(getParam<Real>("rtol")),
    _abs_tol(getParam<Real>("abs_tol")),
    _gtol(getParam<Real>("gtol")),
    _slip_incr_tol(getParam<Real>("slip_incr_tol")),
    _maxiter(getParam<unsigned int>("maxiter")),
    _maxiterg(getParam<unsigned int>("maxitergss")),
    _num_slip_sys_flowrate_props(getParam<unsigned int>("num_slip_sys_flowrate_props")),
    _tan_mod_type(getParam<MooseEnum>("tan_mod_type")),
    _intvar_read_type(getParam<MooseEnum>("intvar_read_type")),
    _num_slip_sys_props(getParam<unsigned int>("num_slip_sys_props")),
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
    _fp(declareADProperty<RankTwoTensor>("fp")), // Plastic deformation gradient
    _fp_old(getMaterialPropertyOld<RankTwoTensor>(
        "fp")), // Plastic deformation gradient of previous increment
    _pk2(declareADProperty<RankTwoTensor>("pk2")), // 2nd Piola Kirchoff Stress
    _pk2_old(getMaterialPropertyOld<RankTwoTensor>(
        "pk2")), // 2nd Piola Kirchoff Stress of previous increment
    _lag_e(declareADProperty<RankTwoTensor>("lage")), // Lagrangian strain
    _lag_e_old(
        getMaterialPropertyOld<RankTwoTensor>("lage")), // Lagrangian strain of previous increment
    _lag_el(declareADProperty<RankTwoTensor>("lagel")), // Lagrangian strain
    _lag_pl(declareADProperty<RankTwoTensor>("lagpl")), // Lagrangian strain
    _gss(declareADProperty<DenseVector<Real>>("gss")),  // Slip system resistances
    _gss_old(getMaterialPropertyOld<DenseVector<Real>>("gss")),
    _acc_slip(declareADProperty<Real>("acc_slip")), // Accumulated slip
    _acc_slip_old(
        getMaterialPropertyOld<Real>("acc_slip")), // Accumulated alip of previous increment
    _update_rot(declareADProperty<RankTwoTensor>(
        "update_rot")), // Rotation tensor considering material rotation and crystal orientation
    _deformation_gradient(getADMaterialProperty<RankTwoTensor>("deformation_gradient")),
    _deformation_gradient_old(getMaterialPropertyOld<RankTwoTensor>("deformation_gradient")),
    _elasticity_tensor_name(_base_name + "elasticity_tensor"),
    _elasticity_tensor(getADMaterialProperty<RankFourTensor>(_elasticity_tensor_name)),
    _crysrot(getADMaterialProperty<RankTwoTensor>("crysrot")),
    _mo(_nss * LIBMESH_DIM),
    _no(_nss * LIBMESH_DIM),
    _a0(_nss),
    _xm(_nss),
    _slip_incr(_nss),
    _tau(_nss),
    _dslipdtau(_nss),
    _s0(_nss),
    _gss_tmp(_nss),
    _gss_tmp_old(_nss),
    //_dgss_dsliprate(_nss, _nss)
    _c(adCoupledValue("c")),
    _c_old(coupledValueOld("c")),
    _kdamage(getParam<Real>("kdamage")),
    _use_current_hist(getParam<bool>("use_current_history_variable")),
    _hist(declareADProperty<Real>("hist")),
    _hist_old(getMaterialPropertyOld<Real>("hist")),
    _Wp(declareADProperty<Real>("Wp")),
    _Wp_old(getMaterialPropertyOld<Real>("Wp")),
    _W0(getParam<Real>("W0")),
    _beta_p(getParam<bool>("beta_p")),
    _beta_e(getParam<bool>("beta_e")),
    _Hall_Petch_const(getParam<Real>("Hall_Petch_const")),
    _grain_size(getParam<Real>("grain_size"))
{
  _err_tol = false;

  if (_num_slip_sys_props > 0)
    _slip_sys_props.resize(_nss * _num_slip_sys_props);

  _pk2_tmp.zero();
  _delta_dfgrd.zero();

  _first_step_iter = false;
  _last_step_iter = false;
  // Initialize variables in the first iteration of substepping
  _first_substep = true;

  _read_from_slip_sys_file = false;
  if (_intvar_read_type == "slip_sys_file")
    _read_from_slip_sys_file = true;

  if (_read_from_slip_sys_file && !(_num_slip_sys_props > 0))
    mooseError("Crystal Plasticity Error: Specify number of internal variable's initial values to "
               "be read from slip system file");

  getSlipSystems();

  if (_slip_sys_flow_prop_file_name.length() != 0)
    readFileFlowRateParams();
  else
    getFlowRateParams();

  if (_slip_sys_hard_prop_file_name.length() != 0)
    readFileHardnessParams();
  else
    getHardnessParams();

  RankTwoTensor::initRandom(_rndm_seed);
}

template <ComputeStage compute_stage>
void
ADFiniteStrainCrystalPlasticityPF<compute_stage>::initQpStatefulProperties()
{
  _stress[_qp].zero();

  _fp[_qp].setToIdentity();

  _pk2[_qp].zero();
  _acc_slip[_qp] = 0.0;
  _lag_e[_qp].zero();

  _update_rot[_qp].setToIdentity();

  _gss[_qp].resize(_nss);
  initSlipSysProps(); // Initializes slip system related properties
  initAdditionalProps();

  _Wp[_qp] = 0;
  _hist[_qp] = 0.0;
}

template <ComputeStage compute_stage>
void
ADFiniteStrainCrystalPlasticityPF<compute_stage>::initSlipSysProps()
{
  switch (_intvar_read_type)
  {
    case 0:
      assignSlipSysRes();
      break;
    case 1:
      readFileInitSlipSysRes();
      break;
    default:
      getInitSlipSysRes();
  }

  // if (_slip_sys_hard_prop_file_name.length() != 0)
  //   readFileHardnessParams();
  // else
  //   getHardnessParams();
}

template <ComputeStage compute_stage>
void
ADFiniteStrainCrystalPlasticityPF<compute_stage>::assignSlipSysRes()
{
  _gss[_qp].resize(_nss);

  for (unsigned int i = 0; i < _nss; ++i)
    _gss[_qp](i) = _slip_sys_props(i) + _Hall_Petch_const / std::sqrt(_grain_size);
}

// Read initial slip system resistances  from .txt file. See test.
template <ComputeStage compute_stage>
void
ADFiniteStrainCrystalPlasticityPF<compute_stage>::readFileInitSlipSysRes()
{
  _gss[_qp].resize(_nss);

  MooseUtils::checkFileReadable(_slip_sys_res_prop_file_name);

  std::ifstream file;
  file.open(_slip_sys_res_prop_file_name.c_str());

  for (unsigned int i = 0; i < _nss; ++i)
  {
    Real v = MetaPhysicL::raw_value(_gss[_qp](i));
    if (!(file >> v))
      mooseError(
          "Error ADFiniteStrainCrystalPlasticityPF: Premature end of slip_sys_res_prop file");
  }

  file.close();
}

// Read initial slip system resistances  from .i file
template <ComputeStage compute_stage>
void
ADFiniteStrainCrystalPlasticityPF<compute_stage>::getInitSlipSysRes()
{
  if (_gprops.size() <= 0)
    mooseError("FiniteStrainCrystalPLasticity: Error in reading slip system resistance properties "
               ": Specify input in .i file or in slip_sys_res_prop_file or in slip_sys_file");

  _gss[_qp].resize(_nss);

  unsigned int num_data_grp = 3; // Number of data per group e.g. start_slip_sys, end_slip_sys,
                                 // value

  for (unsigned int i = 0; i < _gprops.size() / num_data_grp; ++i)
  {
    Real vs, ve;
    unsigned int is, ie;

    vs = _gprops[i * num_data_grp];
    ve = _gprops[i * num_data_grp + 1];

    if (vs <= 0 || ve <= 0)
      mooseError("FiniteStrainCrystalPLasticity: Indices in gss property read must be positive");

    if (vs != floor(vs) || ve != floor(ve))
      mooseError(
          "ADFiniteStrainCrystalPlasticityPF: Error in reading slip system resistances: Values "
          "specifying start and end number of slip system groups should be integer");

    is = static_cast<unsigned int>(vs);
    ie = static_cast<unsigned int>(ve);

    if (is > ie)
      mooseError("FiniteStrainCrystalPLasticity: Start index  should be greater than end index "
                 "in slip system resistance property read");

    for (unsigned int j = is; j <= ie; ++j)
      _gss[_qp](j - 1) = _gprops[i * num_data_grp + 2] + _Hall_Petch_const / std::sqrt(_grain_size);
  }

  for (unsigned int i = 0; i < _nss; ++i)
    if (_gss[_qp](i) <= 0.0)
      mooseError("FiniteStrainCrystalPLasticity: Value of resistance for slip system ",
                 i + 1,
                 " non positive");
}

// Read flow rate parameters from .txt file. See test.
template <ComputeStage compute_stage>
void
ADFiniteStrainCrystalPlasticityPF<compute_stage>::readFileFlowRateParams()
{
  _a0.resize(_nss);
  _xm.resize(_nss);

  MooseUtils::checkFileReadable(_slip_sys_flow_prop_file_name);

  std::ifstream file;
  file.open(_slip_sys_flow_prop_file_name.c_str());

  std::vector<Real> vec;
  vec.resize(_num_slip_sys_flowrate_props);

  for (unsigned int i = 0; i < _nss; ++i)
  {
    for (unsigned int j = 0; j < _num_slip_sys_flowrate_props; ++j)
      if (!(file >> vec[j]))
        mooseError("Error ADFiniteStrainCrystalPlasticityPF: Premature end of "
                   "slip_sys_flow_rate_param file");

    _a0(i) = vec[0];
    _xm(i) = vec[1];
  }

  file.close();
}

// Read flow rate parameters from .i file
template <ComputeStage compute_stage>
void
ADFiniteStrainCrystalPlasticityPF<compute_stage>::getFlowRateParams()
{
  if (_flowprops.size() <= 0)
    mooseError("FiniteStrainCrystalPLasticity: Error in reading flow rate  properties: Specify "
               "input in .i file or a slip_sys_flow_prop_file_name");

  unsigned int num_data_grp = 2 + _num_slip_sys_flowrate_props;

  for (unsigned int i = 0; i < _flowprops.size() / num_data_grp; ++i)
  {
    Real vs, ve;
    unsigned int is, ie;

    vs = _flowprops[i * num_data_grp];
    ve = _flowprops[i * num_data_grp + 1];

    if (vs <= 0 || ve <= 0)
      mooseError("FiniteStrainCrystalPLasticity: Indices in flow rate parameter read must be "
                 "positive integers: is = ",
                 vs,
                 " ie = ",
                 ve);

    if (vs != floor(vs) || ve != floor(ve))
      mooseError("FiniteStrainCrystalPLasticity: Error in reading flow props: Values specifying "
                 "start and end number of slip system groups should be integer");

    is = static_cast<unsigned int>(vs);
    ie = static_cast<unsigned int>(ve);

    if (is > ie)
      mooseError("FiniteStrainCrystalPLasticity: Start index is = ",
                 is,
                 " should be greater than end index ie = ",
                 ie,
                 " in flow rate parameter read");

    for (unsigned int j = is; j <= ie; ++j)
    {
      _a0(j - 1) = _flowprops[i * num_data_grp + 2];
      _xm(j - 1) = _flowprops[i * num_data_grp + 3];
    }
  }

  for (unsigned int i = 0; i < _nss; ++i)
  {
    if (!(_a0(i) > 0.0 && _xm(i) > 0.0))
    {
      mooseWarning("ADFiniteStrainCrystalPlasticityPF: Non-positive flow rate parameters ",
                   _a0(i),
                   ",",
                   _xm(i));
      break;
    }
  }
}

// Read hardness parameters from .txt file
template <ComputeStage compute_stage>
void
ADFiniteStrainCrystalPlasticityPF<compute_stage>::readFileHardnessParams()
{
}

// Read hardness parameters from .i file
template <ComputeStage compute_stage>
void
ADFiniteStrainCrystalPlasticityPF<compute_stage>::getHardnessParams()
{
  if (_hprops.size() <= 0)
    mooseError("FiniteStrainCrystalPLasticity: Error in reading hardness properties: Specify input "
               "in .i file or a slip_sys_hard_prop_file_name");

  _r = _hprops[0];
  _h0 = _hprops[1] + _Hall_Petch_const / std::sqrt(_grain_size);
  _tau_init = _hprops[2];
  _tau_sat = _hprops[3];
}

// Read slip systems from file
template <ComputeStage compute_stage>
void
ADFiniteStrainCrystalPlasticityPF<compute_stage>::getSlipSystems()
{
  Real vec[LIBMESH_DIM];
  std::ifstream fileslipsys;

  MooseUtils::checkFileReadable(_slip_sys_file_name);

  fileslipsys.open(_slip_sys_file_name.c_str());

  for (unsigned int i = 0; i < _nss; ++i)
  {
    // Read the slip normal
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      if (!(fileslipsys >> vec[j]))
        mooseError("Crystal Plasticity Error: Premature end of file reading slip system file \n");

    // Normalize the vectors
    Real mag;
    mag = Utility::pow<2>(vec[0]) + Utility::pow<2>(vec[1]) + Utility::pow<2>(vec[2]);
    mag = std::sqrt(mag);

    for (unsigned j = 0; j < LIBMESH_DIM; ++j)
      _no(i * LIBMESH_DIM + j) = vec[j] / mag;

    // Read the slip direction
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      if (!(fileslipsys >> vec[j]))
        mooseError("Crystal Plasticity Error: Premature end of file reading slip system file\n");

    // Normalize the vectors
    mag = Utility::pow<2>(vec[0]) + Utility::pow<2>(vec[1]) + Utility::pow<2>(vec[2]);
    mag = std::sqrt(mag);

    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      _mo(i * LIBMESH_DIM + j) = vec[j] / mag;

    mag = 0.0;
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      mag += _mo(i * LIBMESH_DIM + j) * _no(i * LIBMESH_DIM + j);

    if (std::abs(mag) > 1e-8)
      mooseError("Crystal Plasicity Error: Slip direction and normal not orthonormal");

    if (_read_from_slip_sys_file)
      for (unsigned int j = 0; j < _num_slip_sys_props; ++j)
        if (!(fileslipsys >> _slip_sys_props(i * _num_slip_sys_props + j)))
          mooseError("Crystal Plasticity Error: Premature end of file reading slip system file");
  }

  fileslipsys.close();
}

// Initialize addtional stateful material properties
template <ComputeStage compute_stage>
void
ADFiniteStrainCrystalPlasticityPF<compute_stage>::initAdditionalProps()
{
}

/**
 * Solves stress residual equation using NR.
 * Updates slip system resistances iteratively.
 */
template <ComputeStage compute_stage>
void
ADFiniteStrainCrystalPlasticityPF<compute_stage>::computeQpStress()
{
  if (isBoundaryMaterial())
    return;

  _gd = Utility::pow<2>(1.0 - _c[_qp]) + _kdamage;
  _gd_old = Utility::pow<2>(1.0 - _c_old[_qp]) + _kdamage;
  //_gd = 1.0;
  //_gd_old = 1.0;

  unsigned int substep_iter = 1; // Depth of substepping; Limited to maximum substep iteration
  unsigned int num_substep = 1;  // Calculated from substep_iter as 2^substep_iter
  Real dt_original = _dt;        // Stores original _dt; Reset at the end of solve
  _first_substep = true;         // Initialize variables at substep_iter = 1

  if (_max_substep_iter > 1)
  {
    _dfgrd_tmp_old = _deformation_gradient_old[_qp];
    if (_dfgrd_tmp_old.det() == 0)
      _dfgrd_tmp_old.addIa(1.0);

    _delta_dfgrd = _deformation_gradient[_qp] - _dfgrd_tmp_old;
    _err_tol = true; // Indicator to continue substepping
  }

  // Substepping loop
  for (substep_iter = 1; substep_iter <= _max_substep_iter; substep_iter++)
  // while (_err_tol && _max_substep_iter > 1)
  {
    _dt = dt_original / num_substep;

    // std::cout << "num_substep = " << num_substep << std::endl;

    for (unsigned int istep = 0; istep < num_substep; ++istep)
    {
      _first_step_iter = false;
      if (istep == 0)
        _first_step_iter = true;

      _last_step_iter = false;
      if (istep == num_substep - 1)
        _last_step_iter = true;

      _dfgrd_scale_factor = (static_cast<Real>(istep) + 1) / num_substep;

      preSolveQp();
      solveQp();

      // if (substep_iter < 4)
      //   _err_tol = true;

      if (_err_tol)
      {
        // substep_iter++;
        // num_substep *= 2;
        // break;
        num_substep *= 2;
        break;
      }
    }

    _first_substep = false; // Prevents reinitialization
    _dt = dt_original;      // Resets dt

#ifdef DEBUG
    if (substep_iter > _max_substep_iter)
      mooseWarning("ADFiniteStrainCrystalPlasticityPF: Failure with substepping");
#endif

    if (!_err_tol || substep_iter > _max_substep_iter)
    {
      postSolveQp(); // Evaluate variables after successful solve or indicate failure
      break;
    }
  }

  // No substepping
  if (_max_substep_iter == 1)
  {
    preSolveQp();
    solveQp();
    postSolveQp();
  }
}

template <ComputeStage compute_stage>
void
ADFiniteStrainCrystalPlasticityPF<compute_stage>::preSolveQp()
{
  // Initialize variable
  if (_first_substep)
  {
    calc_schmid_tensor();
  }

  if (_max_substep_iter == 1)
    _dfgrd_tmp = _deformation_gradient[_qp]; // Without substepping
  else
    _dfgrd_tmp = _dfgrd_scale_factor * _delta_dfgrd + _dfgrd_tmp_old;

  //_dfgrd_tmp.print();

  _err_tol = false;
}

template <ComputeStage compute_stage>
void
ADFiniteStrainCrystalPlasticityPF<compute_stage>::solveQp()
{
  preSolveStatevar();
  solveStatevar();
  if (_err_tol)
    return;
  postSolveStatevar();
}

template <ComputeStage compute_stage>
void
ADFiniteStrainCrystalPlasticityPF<compute_stage>::postSolveQp()
{
  if (_err_tol)
  {
    _err_tol = false;
    if (_gen_rndm_stress_flag)
    {
      if (!_input_rndm_scale_var)
        _rndm_scale_var = MetaPhysicL::raw_value(_elasticity_tensor[_qp](0, 0, 0, 0));

      _stress[_qp] = RankTwoTensor::genRandomSymmTensor(_rndm_scale_var, 1.0);
    }
    // else if (!_fe_problem.getNonlinearSystemBase().computingInitialResidual())
    if (_fe_problem.getNonlinearSystemBase().getCurrentNonlinearIterationNumber() > 0)
      throw MooseException("ADFiniteStrainCrystalPlasticityPF: Constitutive failure");
  }
  else
  {
    // _stress[_qp] = _fe * _pk2[_qp] * _fe.transpose() / _fe.det();
    ADRankTwoTensor iden, ce, ee;
    iden.zero();
    iden.addIa(1.0);
    ce = _fe.transpose() * _fe;
    ee = ce - iden;
    ee *= 0.5;

    ADRankTwoTensor pk2_undamage = _elasticity_tensor[_qp] * ee;

    ADRankTwoTensor D, Q, D_pos, D_neg;
    pk2_undamage.diagonalize(Q, D);

    ADRankTwoTensor pk2_pos, pk2_neg;

    for (unsigned int i = 0; i < 3; ++i)
      D_pos(i, i) = D(i, i) > 0.0 ? D(i, i) : 0;

    pk2_pos = Q * D_pos * Q.transpose();

    ADRankTwoTensor pk2 = _gd * pk2_pos + (_pk2[_qp] - pk2_pos);

    _stress[_qp] = _fe * pk2 * _fe.transpose() / _fe.det();

    // ADRankTwoTensor iden;
    // iden.zero();
    // iden.addIa(1.0);

    _lag_e[_qp] = _deformation_gradient[_qp].transpose() * _deformation_gradient[_qp] - iden;
    _lag_e[_qp] = _lag_e[_qp] * 0.5;

    // ADRankTwoTensor rot;
    // rot = get_current_rotation(_deformation_gradient[_qp]); // Calculate material rotation
    // _update_rot[_qp] = rot * _crysrot[_qp];
  }

  ADRankTwoTensor iden, ce, ee;
  iden.zero();
  iden.addIa(1.0);
  ce = _fe.transpose() * _fe;
  ee = ce - iden;
  ee *= 0.5;

  _lag_el[_qp] = _fe.transpose() * _fe - iden;
  _lag_el[_qp] = _lag_el[_qp] * 0.5;
  _lag_pl[_qp] = _lag_e[_qp] - _lag_el[_qp];

  ADRankTwoTensor pk2_undamage = _elasticity_tensor[_qp] * ee;

  ADRankTwoTensor D, Q, D_pos, D_neg;
  pk2_undamage.diagonalize(Q, D);

  ADRankTwoTensor pk2_pos, pk2_neg;

  for (unsigned int i = 0; i < 3; ++i)
    D_pos(i, i) = D(i, i) > 0.0 ? D(i, i) : 0;

  pk2_pos = Q * D_pos * Q.transpose();

  ADReal G0_pos = 0.0;

  if (_beta_e)
    G0_pos = pk2_pos.doubleContraction(ee) / 2.0;

  DenseVector<ADReal> tau(_nss);

  ADReal plastic_work = 0.0;

  for (unsigned int i = 0; i < _nss; ++i)
    // tau(i) = _pk2[_qp].doubleContraction(_s0[i]);
    tau(i) = pk2_pos.doubleContraction(_s0[i]);

  for (unsigned int i = 0; i < _nss; ++i)
  {
    ADReal gamma = 0.0;
    if (std::abs(_tau(i) / (_gss[_qp](i))) < 1.0e-8)
      gamma = 0.0;
    else if (_tau(i) / (_gss[_qp](i)) > 1.0e-8)
      gamma = _a0(i) * std::pow(std::abs(tau(i) / (_gss[_qp](i))), 1.0 / _xm(i));
    else if (_tau(i) / (_gss[_qp](i)) < 1.0e-8)
      gamma = -_a0(i) * std::pow(std::abs(tau(i) / (_gss[_qp](i))), 1.0 / _xm(i));

    plastic_work += std::abs(gamma * _dt * tau(i));
  }

  _Wp[_qp] = _Wp_old[_qp] + plastic_work;

  if (_beta_p && _Wp[_qp] >= _W0)
    G0_pos += 0.1 * (_Wp[_qp] - _W0);

  // Assign history variable and derivative
  // if (G0_pos > _hist_old[_qp])
  //   _hist[_qp] = G0_pos;
  // else
  //   _hist[_qp] = _hist_old[_qp];

  _hist[_qp] = G0_pos;
}

template <ComputeStage compute_stage>
void
ADFiniteStrainCrystalPlasticityPF<compute_stage>::preSolveStatevar()
{
  if (_max_substep_iter == 1) // No substepping
  {
    _gss_tmp = _gss_old[_qp];
    _accslip_tmp_old = _acc_slip_old[_qp];
  }
  else
  {
    if (_first_step_iter)
    {
      _gss_tmp = _gss_tmp_old = _gss_old[_qp];
      _accslip_tmp_old = _acc_slip_old[_qp];
    }
    else
      _gss_tmp = _gss_tmp_old;
  }
}

template <ComputeStage compute_stage>
void
ADFiniteStrainCrystalPlasticityPF<compute_stage>::solveStatevar()
{
  ADReal gmax;
  // ADReal gdiff;
  unsigned int iterg;
  DenseVector<ADReal> gss_prev(_nss);

  gmax = 1.1 * _gtol;
  iterg = 0;

  while (gmax > _gtol && iterg < _maxiterg) // Check for slip system resistance update tolerance
  {
    preSolveStress();
    solveStress();
    if (_err_tol)
      return;
    postSolveStress();

    gss_prev = _gss_tmp;

    update_slip_system_resistance(); // Update slip system resistance

    gmax = 0.0;
    // for (unsigned i = 0; i < _nss; ++i)
    // {
    //   gdiff = std::abs(gss_prev(i) - _gss_tmp(i)); // Calculate increment size
    //
    //   if (gdiff > gmax)
    //     gmax = gdiff;
    // }
    iterg++;
  }

  if (iterg == _maxiterg)
  {
#ifdef DEBUG
    mooseWarning("FiniteStrainCrystalPLasticity: Hardness Integration error gmax", gmax, "\n");
#endif
    _err_tol = true;
  }
}

template <ComputeStage compute_stage>
void
ADFiniteStrainCrystalPlasticityPF<compute_stage>::postSolveStatevar()
{
  if (_max_substep_iter == 1) // No substepping
  {
    _gss[_qp] = _gss_tmp;
    _acc_slip[_qp] = _accslip_tmp;
  }
  else
  {
    if (_last_step_iter)
    {
      _gss[_qp] = _gss_tmp;
      _acc_slip[_qp] = _accslip_tmp;
    }
    else
    {
      _gss_tmp_old = _gss_tmp;
      _accslip_tmp_old = _accslip_tmp;
    }
  }
}

template <ComputeStage compute_stage>
void
ADFiniteStrainCrystalPlasticityPF<compute_stage>::preSolveStress()
{
  if (_max_substep_iter == 1) // No substepping
  {
    _pk2_tmp = _pk2_old[_qp];
    _fp_old_inv = _fp_old[_qp].inverse();
    _fp_inv = _fp_old_inv;
    _fp_prev_inv = _fp_inv;
  }
  else
  {
    if (_first_step_iter)
    {
      _pk2_tmp = _pk2_tmp_old = _pk2_old[_qp];
      _fp_old_inv = _fp_old[_qp].inverse();
    }
    else
      _pk2_tmp = _pk2_tmp_old;

    _fp_inv = _fp_old_inv;
    _fp_prev_inv = _fp_inv;
  }
}

template <ComputeStage compute_stage>
void
ADFiniteStrainCrystalPlasticityPF<compute_stage>::solveStress()
{
  unsigned int iter = 0;
  ADRankTwoTensor resid, dpk2;
  ADRankFourTensor jac;
  Real rnorm, rnorm0;
  // Real rnorm_prev;

  calc_resid_jacob(resid, jac); // Calculate stress residual
  if (_err_tol)
  {
    //#ifdef DEBUG
    mooseWarning(
        "FiniteStrainCrystalPLasticity: Slip increment exceeds tolerance - Element number ",
        _current_elem->id(),
        " Gauss point = ",
        _qp);
    //#endif
    return;
  }

  rnorm = MetaPhysicL::raw_value(resid.L2norm());
  rnorm0 = rnorm;

  for (iter = 0; iter < _maxiter; ++iter)
  {
    // if (_current_elem->id() == 0 && _qp == 0)
    //   std::cout << "iter = " << iter << std::endl;
    dpk2 = -jac.invSymm() * resid;
    _pk2_tmp = _pk2_tmp + dpk2;
    calc_resid_jacob(resid, jac);
    internalVariableUpdateNRiteration(); // update _fp_prev_inv
    rnorm = MetaPhysicL::raw_value(resid.L2norm());
    if (rnorm < _rtol * rnorm0 || rnorm < _abs_tol)
      break;
  }

  if (iter >= _maxiter)
  {
    //#ifdef DEBUG
    mooseWarning("FiniteStrainCrystalPLasticity: Stress Integration error rmax = ",
                 rnorm,
                 ", rnorm0 = ",
                 rnorm0);
    //#endif
    _err_tol = true;
  }
}

template <ComputeStage compute_stage>
void
ADFiniteStrainCrystalPlasticityPF<compute_stage>::postSolveStress()
{
  if (_max_substep_iter == 1) // No substepping
  {
    _fp[_qp] = _fp_inv.inverse();
    _pk2[_qp] = _pk2_tmp;
  }
  else
  {
    if (_last_step_iter)
    {
      _fp[_qp] = _fp_inv.inverse();
      _pk2[_qp] = _pk2_tmp;
    }
    else
    {
      _fp_old_inv = _fp_inv;
      _pk2_tmp_old = _pk2_tmp;
    }
  }
}

// Update slip system resistance. Overide to incorporate new slip system resistance laws
template <ComputeStage compute_stage>
void
ADFiniteStrainCrystalPlasticityPF<compute_stage>::update_slip_system_resistance()
{
  updateGss();
}

/**
 * Old function to update slip system resistances.
 * Kept to avoid code break at computeQpstress
 */
template <ComputeStage compute_stage>
void
ADFiniteStrainCrystalPlasticityPF<compute_stage>::updateGss()
{
  DenseVector<ADReal> hb(_nss);
  Real qab;

  // DEBUG Real a = _hprops[4]; // Kalidindi
  Real a = _hprops[4];
  // Real a = 2.5;

  _accslip_tmp = _accslip_tmp_old;
  for (unsigned int i = 0; i < _nss; ++i)
    _accslip_tmp += std::abs(_slip_incr(i));

  // Real val = std::cosh(_h0 * _accslip_tmp / (_tau_sat - _tau_init)); // Karthik
  // val = _h0 * std::pow(1.0/val,2.0); // Kalidindi

  for (unsigned int i = 0; i < _nss; ++i)
    // hb(i)=val;
    if ((1.0 - _gss_tmp(i) / _tau_sat) > 1.0e-8)
      hb(i) = _h0 * std::pow(std::abs(1.0 - _gss_tmp(i) / _tau_sat), a);
    else if ((1.0 - _gss_tmp(i) / _tau_sat) < -1.0e-8)
      hb(i) = -_h0 * std::pow(std::abs(1.0 - _gss_tmp(i) / _tau_sat), a);
    else
      hb(i) = 0.0;

  for (unsigned int i = 0; i < _nss; ++i)
  {
    if (_max_substep_iter == 1) // No substepping
      _gss_tmp(i) = _gss_old[_qp](i);
    else
      _gss_tmp(i) = _gss_tmp_old(i);

    for (unsigned int j = 0; j < _nss; ++j)
    {
      unsigned int iplane, jplane;
      iplane = i / 3;
      jplane = j / 3;

      if (iplane == jplane) // Kalidindi
        qab = 1.0;
      else
        qab = _r;

      _gss_tmp(i) += qab * hb(j) * std::abs(_slip_incr(j));
      // _dgss_dsliprate(i, j) =
      //     qab * hb(j) * std::copysign(1.0, MetaPhysicL::raw_value(_slip_incr(j))) * _dt;
    }
  }
}

// Calculates stress residual equation and jacobian
template <ComputeStage compute_stage>
void
ADFiniteStrainCrystalPlasticityPF<compute_stage>::calc_resid_jacob(ADRankTwoTensor & resid,
                                                                   ADRankFourTensor & jac)
{
  calcResidual(resid);
  if (_err_tol)
    return;
  calcJacobian(jac);
}

template <ComputeStage compute_stage>
void
ADFiniteStrainCrystalPlasticityPF<compute_stage>::calcResidual(ADRankTwoTensor & resid)
{
  ADRankTwoTensor iden, ce, ee, ce_pk2, eqv_slip_incr, pk2_new;
  iden.zero();
  iden.addIa(1.0);

  _fe = _dfgrd_tmp * _fp_prev_inv; // _fp_inv  ==> _fp_prev_inv

  ce = _fe.transpose() * _fe;
  ce_pk2 = ce * _pk2_tmp;
  ce_pk2 = ce_pk2 / _fe.det();

  // Calculate Schmid tensor and resolved shear stresses
  for (unsigned int i = 0; i < _nss; ++i)
    _tau(i) = ce_pk2.doubleContraction(_s0[i]);

  getSlipIncrements(); // Calculate dslip,dslipdtau

  if (_err_tol)
    return;

  eqv_slip_incr.zero();
  for (unsigned int i = 0; i < _nss; ++i)
    eqv_slip_incr += _s0[i] * _slip_incr(i);

  eqv_slip_incr = iden - eqv_slip_incr;
  _fp_inv = _fp_old_inv * eqv_slip_incr;
  _fe = _dfgrd_tmp * _fp_inv;

  ce = _fe.transpose() * _fe;
  ee = ce - iden;
  ee *= 0.5;

  // ADRankTwoTensor pk2_undamage = _elasticity_tensor[_qp] * ee;
  // std::vector<ADReal> eigval;
  // ADRankTwoTensor eigvec;
  //
  // pk2_undamage.symmetricEigenvaluesEigenvectors(eigval, eigvec);
  //
  // // Calculate tensors of outerproduct of eigen vectors
  // std::vector<ADRankTwoTensor> etens(LIBMESH_DIM);
  //
  // for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  //   etens[i].vectorOuterProduct(eigvec.column(i), eigvec.column(i));
  //
  // // Separate out positive and negative eigen values
  // std::vector<ADReal> epos(LIBMESH_DIM), eneg(LIBMESH_DIM);
  // for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  // {
  //   epos[i] = (std::abs(eigval[i]) + eigval[i]) / 2.0;
  //   eneg[i] = -(std::abs(eigval[i]) - eigval[i]) / 2.0;
  // }
  //
  // // Calculate the tensile (postive) and compressive (negative) parts of stress
  // ADRankTwoTensor pk2_pos, pk2_neg;
  // pk2_pos.zero();
  // pk2_neg.zero();
  // for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  // {
  //   pk2_pos += etens[i] * epos[i];
  //   pk2_neg += etens[i] * eneg[i];
  // }
  //
  // pk2_new = _gd * pk2_pos + pk2_neg;

  // ADRankTwoTensor pk2_undamage = _elasticity_tensor[_qp] * ee;
  //
  // ADRankTwoTensor D, Q, D_pos, D_neg;
  // pk2_undamage.diagonalize(Q, D);
  //
  // ADRankTwoTensor pk2_pos, pk2_neg;
  //
  // for (unsigned int i = 0; i < 3; ++i)
  //   D_pos(i, i) = D(i, i) > 0.0 ? D(i, i) : 0;
  //
  // D_neg = D - D_pos;
  //
  // pk2_pos = Q * D_pos * Q.transpose();
  // pk2_neg = Q * D_neg * Q.transpose();
  //
  // pk2_new = _gd * pk2_pos + pk2_neg;

  pk2_new = _elasticity_tensor[_qp] * ee;

  resid = _pk2_tmp - pk2_new;
}

template <ComputeStage compute_stage>
void
ADFiniteStrainCrystalPlasticityPF<compute_stage>::calcJacobian(ADRankFourTensor & jac)
{
  ADRankFourTensor dfedfpinv, deedfe, dfpinvdpk2;

  std::vector<ADRankTwoTensor> dtaudpk2(_nss), dfpinvdslip(_nss);

  for (unsigned int i = 0; i < _nss; ++i)
  {
    dtaudpk2[i] = _s0[i];
    dfpinvdslip[i] = -_fp_old_inv * _s0[i];
  }

  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
        dfedfpinv(i, j, k, j) = _dfgrd_tmp(i, k);

  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
      {
        deedfe(i, j, k, i) = deedfe(i, j, k, i) + _fe(k, j) * 0.5;
        deedfe(i, j, k, j) = deedfe(i, j, k, j) + _fe(k, i) * 0.5;
      }

  for (unsigned int i = 0; i < _nss; ++i)
    dfpinvdpk2 += (dfpinvdslip[i] * _dslipdtau(i)).outerProduct(dtaudpk2[i]);

  jac =
      RankFourTensor::IdentityFour() - (_elasticity_tensor[_qp] * deedfe * dfedfpinv * dfpinvdpk2);

  // ADRankTwoTensor ce, ee, iden;
  // iden.zero();
  // iden.addIa(1.0);
  // ce = _fe.transpose() * _fe;
  // ee = ce - iden;
  // ee *= 0.5;
  // ADRankTwoTensor pk2_undamage = _elasticity_tensor[_qp] * ee;
  // ADRankFourTensor Ppos = pk2_undamage.positiveProjectionEigenDecomposition();
  //
  // RankFourTensor I4sym(RankFourTensor::initIdentitySymmetricFour);
  // jac = RankFourTensor::IdentityFour() -
  //       ((_gd * Ppos * _elasticity_tensor[_qp] + (I4sym - Ppos) * _elasticity_tensor[_qp]) *
  //        deedfe * dfedfpinv * dfpinvdpk2);
}

// Calculate slip increment,dslipdtau. Override to modify.
template <ComputeStage compute_stage>
void
ADFiniteStrainCrystalPlasticityPF<compute_stage>::getSlipIncrements()
{
  for (unsigned int i = 0; i < _nss; ++i)
  {
    if (std::abs(_tau(i) / _gss_tmp(i) / 1.0) < 1.0e-8)
      _slip_incr(i) = 0.0;
    else if (_tau(i) / _gss_tmp(i) / 1.0 > 1.0e-8)
      _slip_incr(i) = _a0(i) * std::pow(std::abs(_tau(i) / _gss_tmp(i) / 1.0), 1.0 / _xm(i)) * _dt;
    else if (_tau(i) / _gss_tmp(i) / 1.0 < -1.0e-8)
      _slip_incr(i) = -_a0(i) * std::pow(std::abs(_tau(i) / _gss_tmp(i) / 1.0), 1.0 / _xm(i)) * _dt;

    // std::cout << "_slip_incr(i) = " << _slip_incr(i) << std::endl;

    if (std::abs(_slip_incr(i)) > _slip_incr_tol)
    {
      _err_tol = true;
      //#ifdef DEBUG
      mooseWarning("Maximum allowable slip increment exceeded ", std::abs(_slip_incr(i)));
      //#endif
      return;
    }
  }

  for (unsigned int i = 0; i < _nss; ++i)
    if (std::abs(_tau(i) / _gss_tmp(i) / 1.0) < 1.0e-8)
      _dslipdtau(i) = 0.0;
    else
      _dslipdtau(i) = _a0(i) / _xm(i) *
                      std::pow(std::abs(_tau(i) / _gss_tmp(i) / 1.0), 1.0 / _xm(i) - 1.0) /
                      _gss_tmp(i) / 1.0 * _dt;
}

// // Calls getMatRot to perform RU factorization of a tensor.
// template <ComputeStage compute_stage>
// ADRankTwoTensor
// ADFiniteStrainCrystalPlasticityPF<compute_stage>::get_current_rotation(const ADRankTwoTensor & a)
// {
//   return getMatRot(a);
// }
//
// // Performs RU factorization of a tensor
// template <ComputeStage compute_stage>
// ADRankTwoTensor
// ADFiniteStrainCrystalPlasticityPF<compute_stage>::getMatRot(const ADRankTwoTensor & a)
// {
//   ADRankTwoTensor rot;
//   ADRankTwoTensor c, diag, evec;
//   PetscScalar cmat[LIBMESH_DIM][LIBMESH_DIM], work[10];
//   PetscReal w[LIBMESH_DIM];
//   PetscBLASInt nd = LIBMESH_DIM, lwork = 10, info;
//
//   c = a.transpose() * a;
//
//   for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
//     for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
//       cmat[i][j] = c(i, j);
//
//   LAPACKsyev_("V", "U", &nd, &cmat[0][0], &nd, w, work, &lwork, &info);
//
//   if (info != 0)
//     mooseError("FiniteStrainCrystalPLasticity: DSYEV function call in getMatRot function
//     failed");
//
//   diag.zero();
//
//   for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
//     diag(i, i) = std::sqrt(w[i]);
//
//   for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
//     for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
//       evec(i, j) = cmat[i][j];
//
//   rot = a * ((evec.transpose() * diag * evec).inverse());
//
//   return rot;
// }

template <ComputeStage compute_stage>
void
ADFiniteStrainCrystalPlasticityPF<compute_stage>::calc_schmid_tensor()
{
  DenseVector<Real> mo(LIBMESH_DIM * _nss), no(LIBMESH_DIM * _nss);

  // Update slip direction and normal with crystal orientation
  for (unsigned int i = 0; i < _nss; ++i)
  {
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
    {
      mo(i * LIBMESH_DIM + j) = 0.0;
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
        mo(i * LIBMESH_DIM + j) =
            mo(i * LIBMESH_DIM + j) + _crysrot[_qp](j, k) * _mo(i * LIBMESH_DIM + k);
    }

    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
    {
      no(i * LIBMESH_DIM + j) = 0.0;
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
        no(i * LIBMESH_DIM + j) =
            no(i * LIBMESH_DIM + j) + _crysrot[_qp](j, k) * _no(i * LIBMESH_DIM + k);
    }
  }

  // Calculate Schmid tensor and resolved shear stresses
  for (unsigned int i = 0; i < _nss; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
        _s0[i](j, k) = mo(i * LIBMESH_DIM + j) * no(i * LIBMESH_DIM + k);
}

template <ComputeStage compute_stage>
bool
ADFiniteStrainCrystalPlasticityPF<compute_stage>::line_search_update(const ADReal rnorm_prev,
                                                                     const ADRankTwoTensor dpk2)
{
  if (_lsrch_method == "CUT_HALF")
  {
    ADReal rnorm;
    ADRankTwoTensor resid;
    Real step = 1.0;

    do
    {
      _pk2_tmp = _pk2_tmp - step * dpk2;
      step /= 2.0;
      _pk2_tmp = _pk2_tmp + step * dpk2;

      calcResidual(resid);
      rnorm = resid.L2norm();
    } while (rnorm > rnorm_prev && step > _min_lsrch_step);

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
    ADReal s_m = 1000.0;
    ADReal rnorm = 1000.0;

    ADRankTwoTensor resid;
    calcResidual(resid);
    ADReal s_b = resid.doubleContraction(dpk2);
    ADReal rnorm1 = resid.L2norm();
    _pk2_tmp = _pk2_tmp - dpk2;
    calcResidual(resid);
    ADReal s_a = resid.doubleContraction(dpk2);
    ADReal rnorm0 = resid.L2norm();
    _pk2_tmp = _pk2_tmp + dpk2;

    if ((rnorm1 / rnorm0) < _lsrch_tol || s_a * s_b > 0)
    {
      calcResidual(resid);
      return true;
    }

    while ((rnorm / rnorm0) > _lsrch_tol && count < _lsrch_max_iter)
    {
      _pk2_tmp = _pk2_tmp - step * dpk2;
      step = 0.5 * (step_b + step_a);
      _pk2_tmp = _pk2_tmp + step * dpk2;
      calcResidual(resid);
      s_m = resid.doubleContraction(dpk2);
      rnorm = resid.L2norm();

      if (s_m * s_a < 0.0)
      {
        step_b = step;
        s_b = s_m;
      }
      if (s_m * s_b < 0.0)
      {
        step_a = step;
        s_a = s_m;
      }
      count++;
    }

    if ((rnorm / rnorm0) < _lsrch_tol && count < _lsrch_max_iter)
      return true;

    return false;
  }
  else
  {
    mooseError("Line search meothod is not provided.");
    return false;
  }
}

template <ComputeStage compute_stage>
void
ADFiniteStrainCrystalPlasticityPF<compute_stage>::internalVariableUpdateNRiteration()
{
  _fp_prev_inv = _fp_inv; // update _fp_prev_inv
}
