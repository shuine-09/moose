//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeLinearElasticPFFractureStress.h"

registerMooseObject("TensorMechanicsApp", ComputeLinearElasticPFFractureStress);

template <>
InputParameters
validParams<ComputeLinearElasticPFFractureStress>()
{
  InputParameters params = validParams<ComputeIsotropicLinearElasticPFFractureStress>();
  params.addClassDescription("Computes the stress and free energy derivatives for the phase field "
                             "fracture model, with linear anistropic elasticity");
  return params;
}

ComputeLinearElasticPFFractureStress::ComputeLinearElasticPFFractureStress(
    const InputParameters & parameters)
  : ComputeIsotropicLinearElasticPFFractureStress(parameters)
{
}

void
ComputeLinearElasticPFFractureStress::computeQpStress()
{
  Real c = _c[_qp];

  // Compute eigenvectors and eigenvalues of stress
  RankTwoTensor eigvec;
  std::vector<Real> eigval(LIBMESH_DIM);
  RankTwoTensor stress = _elasticity_tensor[_qp] * _mechanical_strain[_qp];
  RankFourTensor I4sym(RankFourTensor::initIdentitySymmetricFour);

  RankTwoTensor stress0pos, stress0neg;

  // if (c < 0.9)
  {
    // projection tensor
    RankFourTensor proj_pos = stress.positiveProjectionEigenDecomposition(eigval, eigvec);
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

    // Calculate the tensile (postive) and compressive (negative) parts of stress
    for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    {
      stress0pos += etens[i] * epos[i];
      stress0neg += etens[i] * eneg[i];
    }

    // Energy with positive principal stress
    const Real G0_pos = (stress0pos).doubleContraction(_mechanical_strain[_qp]) / 2.0;
    const Real G0_neg = (stress0neg).doubleContraction(_mechanical_strain[_qp]) / 2.0;

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

    _Jacobian_mult[_qp] = ((1.0 - c) * (1.0 - c) * (1.0 - _kdamage) + _kdamage) *
                              (proj_pos * _elasticity_tensor[_qp]) +
                          proj_neg * _elasticity_tensor[_qp];
  }
  /*
  else
  {
    RankTwoTensor I(RankTwoTensor::initIdentity);

    if (stress.trace() > 0)
      stress0pos = 1.0 / LIBMESH_DIM * stress.trace() * I + stress.deviatoric();
    else
      stress0pos = stress.deviatoric();

    stress0neg = stress - stress0pos;

    // Energy with positive principal stress
    const Real G0_pos = (stress0pos).doubleContraction(_mechanical_strain[_qp]) / 2.0;
    const Real G0_neg = (stress0neg).doubleContraction(_mechanical_strain[_qp]) / 2.0;

    // Assign history variable and derivative
    if (G0_pos > _hist_old[_qp])
      _hist[_qp] = G0_pos;
    else
      _hist[_qp] = _hist_old[_qp];

    Real hist_variable = _hist_old[_qp];
    if (_use_current_hist)
      hist_variable = _hist[_qp];

    _stress[_qp] = stress0pos * ((1.0 - c) * (1.0 - c) * (1 - _kdamage) + _kdamage) + stress0neg;

    // Elastic free energy density
    _F[_qp] = hist_variable * ((1.0 - c) * (1.0 - c) * (1 - _kdamage) + _kdamage) + G0_neg +
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

    Real gd = (1.0 - c) * (1.0 - c) * (1.0 - _kdamage) + _kdamage;

    if (stress.trace() > 0)
      _Jacobian_mult[_qp] = gd * _elasticity_tensor[_qp];
    else
      _Jacobian_mult[_qp] = (gd * (I4sym - 1.0 / LIBMESH_DIM * I.outerProduct(I)) +
                             1.0 / LIBMESH_DIM * I.outerProduct(I)) *
                            _elasticity_tensor[_qp];

    RankTwoTensor test = _stress[_qp] - _Jacobian_mult[_qp] * _mechanical_strain[_qp];
    //if (test.L2norm() > 1e-8)
    //  mooseWarning("Wrong Jacobian");
  }
  */
}
