//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeIsotropicLinearElasticPFFractureStress.h"
#include "MathUtils.h"

registerMooseObject("TensorMechanicsApp", ComputeIsotropicLinearElasticPFFractureStress);

template <>
InputParameters
validParams<ComputeIsotropicLinearElasticPFFractureStress>()
{
  InputParameters params = validParams<ComputeStressBase>();
  params.addClassDescription("Computes the stress and free energy derivatives for the phase field "
                             "fracture model, with linear isotropic elasticity");
  params.addRequiredCoupledVar("c", "Order parameter for damage");
  params.addParam<Real>("kdamage", 0.0, "Stiffness of damaged matrix");
  params.addParam<bool>(
      "use_current_history_variable", false, "Use the current value of the history variable.");
  params.addParam<MaterialPropertyName>(
      "F_name", "E_el", "Name of material property storing the elastic energy");
  params.addParam<bool>(
      "use_spectral",
      true,
      "True: use spectrual decomposition; False: use volumetric/deviatoric decompostion");
  params.addParam<bool>("use_vi_solver", false, "Use vi solver.");
  return params;
}

ComputeIsotropicLinearElasticPFFractureStress::ComputeIsotropicLinearElasticPFFractureStress(
    const InputParameters & parameters)
  : ComputeStressBase(parameters),
    _c(coupledValue("c")),
    _kdamage(getParam<Real>("kdamage")),
    _use_current_hist(getParam<bool>("use_current_history_variable")),
    _l(getMaterialProperty<Real>("l")),
    _gc(getMaterialProperty<Real>("gc_prop")),
    _F(declareProperty<Real>(getParam<MaterialPropertyName>("F_name"))),
    _dFdc(declarePropertyDerivative<Real>(getParam<MaterialPropertyName>("F_name"),
                                          getVar("c", 0)->name())),
    _d2Fdc2(declarePropertyDerivative<Real>(
        getParam<MaterialPropertyName>("F_name"), getVar("c", 0)->name(), getVar("c", 0)->name())),
    _dstress_dc(
        declarePropertyDerivative<RankTwoTensor>(_base_name + "stress", getVar("c", 0)->name())),
    _d2Fdcdstrain(declareProperty<RankTwoTensor>("d2Fdcdstrain")),
    _hist(declareProperty<Real>("hist")),
    _hist_old(getMaterialPropertyOld<Real>("hist")),
    _use_spectral(getParam<bool>("use_spectral")),
    _use_vi_solver(getParam<bool>("use_vi_solver"))
{
}

void
ComputeIsotropicLinearElasticPFFractureStress::initQpStatefulProperties()
{
  _hist[_qp] = 0.0;
}

void
ComputeIsotropicLinearElasticPFFractureStress::computeQpStress()
{
  if (_use_spectral)
  {
    const Real c = _c[_qp];

    // Isotropic elasticity is assumed and should be enforced
    const Real lambda = _elasticity_tensor[_qp](0, 0, 1, 1);
    const Real mu = _elasticity_tensor[_qp](0, 1, 0, 1);

    // Compute eigenvectors and eigenvalues of mechanical strain and projection tensor
    RankTwoTensor eigvec;
    std::vector<Real> eigval(LIBMESH_DIM);
    RankFourTensor proj_pos =
        _mechanical_strain[_qp].positiveProjectionEigenDecomposition(eigval, eigvec);
    RankFourTensor I4sym(RankFourTensor::initIdentitySymmetricFour);
    RankFourTensor proj_neg = I4sym - proj_pos;

    // Calculate tensors of outerproduct of eigen vectors
    std::vector<RankTwoTensor> etens(LIBMESH_DIM);

    for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
      for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
        for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
          etens[i](j, k) = eigvec(j, i) * eigvec(k, i);

    // Separate out positive and negative eigen values
    std::vector<Real> epos(LIBMESH_DIM), eneg(LIBMESH_DIM);
    for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    {
      epos[i] = (std::abs(eigval[i]) + eigval[i]) / 2.0;
      eneg[i] = (std::abs(eigval[i]) - eigval[i]) / 2.0;
    }

    // Seprate positive and negative sums of all eigenvalues
    Real etr = 0.0;
    for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
      etr += eigval[i];

    const Real etrpos = (std::abs(etr) + etr) / 2.0;
    const Real etrneg = (std::abs(etr) - etr) / 2.0;

    // Calculate the tensile (postive) and compressive (negative) parts of stress
    RankTwoTensor stress0pos, stress0neg;
    for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    {
      stress0pos += etens[i] * (lambda * etrpos + 2.0 * mu * epos[i]);
      stress0neg += etens[i] * (lambda * etrneg + 2.0 * mu * eneg[i]);
    }

    // sum squares of epos and eneg
    Real pval(0.0), nval(0.0);
    for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    {
      pval += epos[i] * epos[i];
      nval += eneg[i] * eneg[i];
    }

    // Energy with positive principal strains
    const Real G0_pos = lambda * etrpos * etrpos / 2.0 + mu * pval;
    const Real G0_neg = lambda * etrneg * etrneg / 2.0 + mu * nval;

    // Assign history variable and derivative
    // Assign history variable and derivative
    if (_use_vi_solver)
    {
      _hist[_qp] = G0_pos;
    }
    else
    {
      if (G0_pos > _hist_old[_qp])
        _hist[_qp] = G0_pos;
      else
        _hist[_qp] = _hist_old[_qp];
    }

    Real hist_variable = _hist_old[_qp];
    if (_use_current_hist)
      hist_variable = _hist[_qp];

    // Damage associated with positive component of stress
    _stress[_qp] = stress0pos * ((1.0 - c) * (1.0 - c) * (1 - _kdamage) + _kdamage) - stress0neg;

    // Elastic free energy density
    _F[_qp] = hist_variable * ((1.0 - c) * (1.0 - c) * (1 - _kdamage) + _kdamage) - G0_neg +
              _gc[_qp] / (2 * _l[_qp]) * c * c;

    // derivative of elastic free energy density wrt c
    _dFdc[_qp] = -hist_variable * 2.0 * (1.0 - c) * (1 - _kdamage) + _gc[_qp] / _l[_qp] * c;

    // 2nd derivative of elastic free energy density wrt c
    _d2Fdc2[_qp] = hist_variable * 2.0 * (1 - _kdamage) + _gc[_qp] / _l[_qp];

    // 2nd derivative wrt c and strain = 0.0 if we used the previous step's history varible
    if (_use_current_hist)
      _d2Fdcdstrain[_qp] = -stress0pos * 2.0 * (1.0 - c) * (1 - _kdamage);

    // Used in StressDivergencePFFracTensors off-diagonal Jacobian
    _dstress_dc[_qp] = -stress0pos * 2.0 * (1.0 - c) * (1 - _kdamage);

    // Calculate dstress/dstrain
    // dstress/dstrain = dstress_pos/dstrain_pos * dstrain_pos/dstrain
    //                 + dstress_neg/dstrain_neg * dstrain_neg/dstrain
    // dstrain_pos/dstrain = proj_pos (fourth order projection tensor)
    // dstrain_neg/dstrain = proj_neg (fourth order projection tensor)
    // Note that proj_pos + proj_neg = symmetric identify fourth order tensor

    // d trace(A) / dA
    RankTwoTensor I(RankTwoTensor::initIdentity);
    RankFourTensor dtraceAdA = I.outerProduct(I);

    _Jacobian_mult[_qp] = ((1.0 - c) * (1.0 - c) * (1 - _kdamage) + _kdamage) *
                              (lambda * dtraceAdA * MathUtils::heavyside(etr) + 2 * mu * proj_pos) +
                          (lambda * dtraceAdA * MathUtils::heavyside(-etr) + 2 * mu * proj_neg);
  }
  else
  {

    const Real c = _c[_qp];

    // Compute cfactor to have c > 1 behave the same as c = 1
    Real cfactor = 0.0;
    if (c < 1.0)
      cfactor = 1.0;

    // Isotropic elasticity is assumed and should be enforced
    const Real lambda = _elasticity_tensor[_qp](0, 0, 1, 1);
    const Real mu = _elasticity_tensor[_qp](0, 1, 0, 1);
    const Real k = lambda + 2.0 * mu / LIBMESH_DIM;

    RankFourTensor I4(RankFourTensor::initIdentitySymmetricFour);
    RankTwoTensor I2(RankTwoTensor::initIdentity);
    RankFourTensor I2I2 = I2.outerProduct(I2);
    RankFourTensor Jacobian_active, Jacobian_passive;
    RankTwoTensor strain0vol, strain0dev;
    RankTwoTensor stress0pos, stress0neg;
    Real strain0tr, strain0tr_neg, strain0tr_pos;

    strain0dev = _mechanical_strain[_qp].deviatoric();
    strain0vol = _mechanical_strain[_qp] - strain0dev;
    strain0tr = _mechanical_strain[_qp].trace();
    strain0tr_neg = std::min(strain0tr, 0.0);
    strain0tr_pos = strain0tr - strain0tr_neg;

    stress0neg = k * strain0tr_neg * I2;
    stress0pos = _elasticity_tensor[_qp] * _mechanical_strain[_qp] - stress0neg;

    // Energy with positive principal strains
    RankTwoTensor strain0dev2 = strain0dev * strain0dev;
    const Real G0_pos = 0.5 * k * strain0tr_pos * strain0tr_pos + mu * strain0dev2.trace();
    const Real G0_neg = 0.5 * k * strain0tr_neg * strain0tr_neg;

    // Assign history variable and derivative
    if (_use_vi_solver)
    {
      _hist[_qp] = G0_pos;
    }
    else
    {
      if (G0_pos > _hist_old[_qp])
        _hist[_qp] = G0_pos;
      else
        _hist[_qp] = _hist_old[_qp];
    }

    Real hist_variable = _hist_old[_qp];
    if (_use_current_hist)
      hist_variable = _hist[_qp];

    // Elastic free energy density
    _F[_qp] = hist_variable * cfactor * ((1.0 - c) * (1.0 - c) * (1 - _kdamage) + _kdamage) -
              G0_neg + _gc[_qp] / (2 * _l[_qp]) * c * c;

    // derivative of elastic free energy density wrt c
    _dFdc[_qp] =
        -hist_variable * 2.0 * cfactor * (1.0 - c) * (1 - _kdamage) + _gc[_qp] / _l[_qp] * c;

    // 2nd derivative of elastic free energy density wrt c
    _d2Fdc2[_qp] = hist_variable * 2.0 * cfactor * (1 - _kdamage) + _gc[_qp] / _l[_qp];

    // 2nd derivative wrt c and strain. Note that I am ignoring the history variable here, but
    // this approach gets the fastest convergence
    _d2Fdcdstrain[_qp] = -stress0pos * cfactor * (1.0 - c) * (1 - _kdamage);

    _stress[_qp] =
        (cfactor * ((1.0 - c) * (1.0 - c) * (1 - _kdamage) + _kdamage)) * stress0pos + stress0neg;

    _dstress_dc[_qp] = -stress0pos * 2.0 * cfactor * (1.0 - c) * (1 - _kdamage);

    // jacobian
    if (strain0tr < 0)
      Jacobian_passive = k * I2I2;
    Jacobian_active = _elasticity_tensor[_qp] - Jacobian_passive;
    _Jacobian_mult[_qp] =
        (cfactor * ((1.0 - c) * (1.0 - c) * (1 - _kdamage) + _kdamage)) * Jacobian_active +
        Jacobian_passive;
  }
}
