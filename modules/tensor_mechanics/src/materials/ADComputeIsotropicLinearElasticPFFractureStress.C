//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADComputeIsotropicLinearElasticPFFractureStress.h"

registerADMooseObject("TensorMechanicsApp", ADComputeIsotropicLinearElasticPFFractureStress);

defineADValidParams(ADComputeIsotropicLinearElasticPFFractureStress,
                    ADComputeStressBase,
                    params.addClassDescription(
                        "Computes the stress and free energy derivatives for the phase field "
                        "fracture model, with linear isotropic elasticity");
                    params.addRequiredCoupledVar("c", "Order parameter for damage");
                    params.addParam<Real>("kdamage", 0.0, "Stiffness of damaged matrix");
                    params.addParam<bool>("use_current_history_variable",
                                          false,
                                          "Use the current value of the history variable."););

template <ComputeStage compute_stage>
ADComputeIsotropicLinearElasticPFFractureStress<compute_stage>::
    ADComputeIsotropicLinearElasticPFFractureStress(const InputParameters & parameters)
  : ADComputeStressBase<compute_stage>(parameters),
    _c(adCoupledValue("c")),
    _c_old(coupledValueOld("c")),
    _kdamage(adGetParam<Real>("kdamage")),
    _use_current_hist(adGetParam<bool>("use_current_history_variable")),
    _hist(adDeclareADProperty<Real>("hist")),
    _hist_old(adGetMaterialPropertyOld<Real>("hist"))
{
}

template <ComputeStage compute_stage>
void
ADComputeIsotropicLinearElasticPFFractureStress<compute_stage>::computeQpStress()
{
  const ADReal c = _c[_qp];

  // Zero out values when c > 1
  ADReal cfactor = 1.0;
  if (c > 1.0)
    cfactor = 0.0;

  // const ADReal lambda = _elasticity_tensor[_qp](0, 0, 1, 1);
  // const ADReal mu = _elasticity_tensor[_qp](0, 1, 0, 1);
  // const ADReal k = lambda + 2.0 * mu / LIBMESH_DIM;
  //
  // ADRankTwoTensor I2;
  // I2.addIa(1.0);
  // ADRankTwoTensor strain0vol, strain0dev;
  // ADRankTwoTensor stress0pos, stress0neg;
  // ADReal strain0tr, strain0tr_neg, strain0tr_pos;
  //
  // strain0dev = _mechanical_strain[_qp].deviatoric();
  // strain0vol = _mechanical_strain[_qp] - strain0dev;
  // strain0tr = _mechanical_strain[_qp].trace();
  // strain0tr_neg = std::min(strain0tr, 0.0);
  // strain0tr_pos = strain0tr - strain0tr_neg;
  //
  // stress0neg = k * strain0tr_neg * I2;
  // stress0pos = _elasticity_tensor[_qp] * _mechanical_strain[_qp] - stress0neg;
  //
  // ADRankTwoTensor strain0dev2 = strain0dev * strain0dev;
  // const ADReal G0_pos = 0.5 * k * strain0tr_pos * strain0tr_pos + mu * strain0dev2.trace();
  // // const ADReal G0_neg = 0.5 * k * strain0tr_neg * strain0tr_neg;
  //
  // // Assign history variable and derivative
  // if (G0_pos > _hist_old[_qp])
  //   _hist[_qp] = G0_pos;
  // else
  //   _hist[_qp] = _hist_old[_qp];
  //
  // _stress[_qp] =
  //     (cfactor * ((1.0 - c) * (1.0 - c) * (1 - _kdamage) + _kdamage)) * stress0pos + stress0neg;

  //===================================================================================================
  // Isotropic elasticity is assumed and should be enforced
  // const ADReal lambda = _elasticity_tensor[_qp](0, 0, 1, 1);
  // const ADReal mu = _elasticity_tensor[_qp](0, 1, 0, 1);
  //
  // std::array<ADReal, 3> eval;
  // std::array<std::array<ADReal, 3>, 3> evec;
  //
  // // Calculate tensors of outerproduct of eigen vectors
  // std::vector<ADRankTwoTensor> etens(LIBMESH_DIM);
  //
  // _mechanical_strain[_qp].calculateEigenvalueAndVector(eval, evec);
  //
  // // ADRankTwoTensor test;
  // // test(0, 0) = 1;
  // // test(0, 1) = 2;
  // // test(0, 2) = 3;
  // // test(1, 0) = 2;
  // // test(1, 1) = 4;
  // // test(1, 2) = 5;
  // // test(2, 1) = 5;
  // // test(2, 2) = 4;
  // // test(2, 2) = 6;
  // //
  // // test.calculate_eigen(eval, evec);
  // //
  // // std::cout << "eigval = " << eval[0] << ", " << eval[1] << ", " << eval[2] << std::endl;
  // //
  // // std::cout << "eigvec 0 = " << evec[0][0] << ", " << evec[0][1] << ", " << evec[0][2] <<
  // // std::endl; std::cout << "eigvec 1 = " << evec[1][0] << ", " << evec[1][1] << ", " << evec[1][2]
  // // << std::endl; std::cout << "eigvec 2 = " << evec[2][0] << ", " << evec[2][1] << ", " <<
  // // evec[2][2] << std::endl;
  // //
  // // mooseError("STOP");
  //
  // for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  //   for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
  //     for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
  //       //      etens[i](j, k) = evec[j][i] * evec[k][i];
  //       etens[i](j, k) = evec[i][j] * evec[i][k];
  //
  // // Separate out positive and negative eigen values
  // std::vector<ADReal> epos(LIBMESH_DIM), eneg(LIBMESH_DIM);
  // for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  // {
  //   epos[i] = (std::abs(eval[i]) + eval[i]) / 2.0;
  //   eneg[i] = (std::abs(eval[i]) - eval[i]) / 2.0;
  // }
  //
  // // Seprate positive and negative sums of all eigenvalues
  // ADReal etr = 0.0;
  // for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  //   etr += eval[i];
  //
  // const ADReal etrpos = (std::abs(etr) + etr) / 2.0;
  // const ADReal etrneg = (std::abs(etr) - etr) / 2.0;
  //
  // // Calculate the tensile (postive) and compressive (negative) parts of stress
  // ADRankTwoTensor stress0pos, stress0neg;
  // for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  // {
  //   stress0pos += etens[i] * (lambda * etrpos + 2.0 * mu * epos[i]);
  //   stress0neg += etens[i] * (lambda * etrneg + 2.0 * mu * eneg[i]);
  // }
  //
  // // sum squares of epos and eneg
  // ADReal pval(0.0), nval(0.0);
  // for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  // {
  //   pval += epos[i] * epos[i];
  //   nval += eneg[i] * eneg[i];
  // }
  //
  // // Energy with positive principal strains
  // const ADReal G0_pos = lambda * etrpos * etrpos / 2.0 + mu * pval;
  // // const ADReal G0_neg = lambda * etrneg * etrneg / 2.0 + mu * nval;
  //
  // // Assign history variable and derivative
  // if (G0_pos > _hist_old[_qp])
  //   _hist[_qp] = G0_pos;
  // else
  //   _hist[_qp] = _hist_old[_qp];
  //
  // // Damage associated with positive component of stress
  // _stress[_qp] = stress0pos * ((1.0 - c) * (1.0 - c) * (1 - _kdamage) + _kdamage) - stress0neg;

  //==========================================================================

  ADRankTwoTensor stress = _elasticity_tensor[_qp] * _mechanical_strain[_qp];

  // std::array<ADReal, 3> eval;
  // std::array<std::array<ADReal, 3>, 3> evec;
  //
  // stress.calculateEigenvalueAndVector(eval, evec);
  //
  // // Project the positive and negative stresses
  // ADRankTwoTensor stress0pos, stress0neg, eigen_value_r2, eigen_vec_r2;
  //
  // for (unsigned int i = 0; i < 3; ++i)
  //   for (unsigned int j = 0; j < 3; ++j)
  //   {
  //     if (i == j)
  //       eigen_value_r2(i, i) = eval[i];
  //     eigen_vec_r2(i, j) = evec[i][j];
  //   }
  //
  // stress0pos = eigen_vec_r2.transpose() * eigen_value_r2 * eigen_vec_r2;

  ADRankTwoTensor D, Q, D_pos, D_neg;

  stress.diagonalize(Q, D);

  ADRankTwoTensor stress0pos, stress0neg;

  for (unsigned int i = 0; i < 3; ++i)
    D_pos(i, i) = D(i, i) > 0.0 ? D(i, i) : 0;

  D_neg = D - D_pos;

  stress0pos = Q * D_pos * Q.transpose();
  stress0neg = Q * D_neg * Q.transpose();

  // stress = Q*D*QT
  // std::cout << "stress = " << stress << std::endl;
  // std::cout << "Q*D*QT = " << Q * D * Q.transpose() << std::endl;

  // Compute the positive and negative elastic energies
  ADReal G0_pos = (stress0pos).doubleContraction(_mechanical_strain[_qp]) / 2.0;

  // Assign history variable and derivative
  if (G0_pos > _hist_old[_qp])
    _hist[_qp] = G0_pos;
  else
    _hist[_qp] = _hist_old[_qp];

  // Damage associated with positive component of stress
  _stress[_qp] = stress0pos * ((1.0 - c) * (1.0 - c) * (1 - _kdamage) + _kdamage) + stress0neg;

  // _stress[_qp] = (1.0 - c) * (1.0 - c) * _elasticity_tensor[_qp] * _mechanical_strain[_qp];
}
