//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeMultiPhaseLinearElasticPFFractureStress.h"

registerMooseObject("TensorMechanicsApp", ComputeMultiPhaseLinearElasticPFFractureStress);

template <>
InputParameters
validParams<ComputeMultiPhaseLinearElasticPFFractureStress>()
{
  InputParameters params = validParams<ComputeLinearElasticPFFractureStress>();
  params.addParam<std::vector<MaterialPropertyName>>(
      "h", "Switching Function Materials that provide h(eta_i)");
  params.addRequiredParam<std::vector<std::string>>("phase_base",
                                                    "Base names for the Phase strains");
  params.addParam<std::string>("base_name", "Base name for the computed global stress (optional)");
  params.addParam<Real>("pressure", 0.0, "Pressure on the crack surfaces.");
  params.addClassDescription("Computes the stress and free energy derivatives for the phase field "
                             "fracture model, with linear anistropic elasticity");
  params.addRequiredCoupledVar("bubble", "bubble");
  return params;
}

ComputeMultiPhaseLinearElasticPFFractureStress::ComputeMultiPhaseLinearElasticPFFractureStress(
    const InputParameters & parameters)
  : ComputeLinearElasticPFFractureStress(parameters),
    _h_list(getParam<std::vector<MaterialPropertyName>>("h")),
    _n_phase(_h_list.size()),
    _h_eta(_n_phase),
    _phase_base(getParam<std::vector<std::string>>("phase_base")),
    _phase_stress(_n_phase),
    _dphase_stress_dstrain(_n_phase),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _pressure(getParam<Real>("pressure")),
    _bubble(coupledValue("bubble"))
{
  // verify parameter length
  if (_n_phase != _phase_base.size())
    mooseError(
        "h and phase_base input vectors need to have the same length in MultiPhaseStressMaterial ",
        name());

  for (unsigned int i = 0; i < _n_phase; ++i)
  {
    _h_eta[i] = &getMaterialProperty<Real>(_h_list[i]);
    _phase_stress[i] = &getMaterialProperty<RankTwoTensor>(_phase_base[i] + "_stress");
    _dphase_stress_dstrain[i] =
        &getMaterialProperty<RankFourTensor>(_phase_base[i] + "_Jacobian_mult");
  }
}

void
ComputeMultiPhaseLinearElasticPFFractureStress::computeQpStress()
{
  Real c = _c[_qp];

  // Compute eigenvectors and eigenvalues of stress
  RankTwoTensor eigvec;
  std::vector<Real> eigval(LIBMESH_DIM);

  RankTwoTensor stress;
  stress.zero();

  stress = _bubble[_qp] * (*_phase_stress[0])[_qp] + (1 - _bubble[_qp]) * (*_phase_stress[1])[_qp];

  // projection tensor
  RankFourTensor proj_pos = stress.positiveProjectionEigenDecomposition(eigval, eigvec);
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

  // Calculate the tensile (postive) and compressive (negative) parts of stress
  RankTwoTensor stress0pos, stress0neg;
  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  {
    stress0pos += etens[i] * epos[i];
    stress0neg += etens[i] * eneg[i];
  }

  // Energy with positive principal stress
  const Real G0_pos = (stress0pos).doubleContraction(_mechanical_strain[_qp]) / 2.0;
  const Real G0_neg = (stress0neg).doubleContraction(_mechanical_strain[_qp]) / 2.0;

  // Assign history variable and derivative
//  if (G0_pos > _hist_old[_qp])
    _hist[_qp] = G0_pos;
//  else
//    _hist[_qp] = _hist_old[_qp];

  Real hist_variable = _hist_old[_qp];
  if (_use_current_hist)
    hist_variable = _hist[_qp];

  _stress[_qp] = stress0pos * ((1.0 - c) * (1.0 - c) * (1 - _kdamage) + _kdamage) - stress0neg +
                 _pressure * RankTwoTensor(RankTwoTensor::initIdentity) * c * c;

  // Elastic free energy density
  _F[_qp] = hist_variable * ((1.0 - c) * (1.0 - c) * (1 - _kdamage) + _kdamage) - G0_neg +
            _gc_prop[_qp] / (2 * _l[_qp]) * c * c +
            c * c * _mechanical_strain[_qp].trace() * _pressure;

  // derivative of elastic free energy density wrt c
  _dFdc[_qp] = -hist_variable * 2.0 * (1.0 - c) * (1 - _kdamage) + _gc_prop[_qp] / _l[_qp] * c +
               2 * c * _mechanical_strain[_qp].trace() * _pressure;

  // 2nd derivative of elastic free energy density wrt c
  _d2Fdc2[_qp] = hist_variable * 2.0 * (1 - _kdamage) + _gc_prop[_qp] / _l[_qp] +
                 2 * _mechanical_strain[_qp].trace() * _pressure;

  // 2nd derivative wrt c and strain = 0.0 if we used the previous step's history varible
  // if (_use_current_hist)
  _d2Fdcdstrain[_qp] = 2 * c * _pressure * RankTwoTensor(RankTwoTensor::initIdentity);

  // Used in StressDivergencePFFracTensors off-diagonal Jacobian
  // _dstress_dc[_qp] = -stress0pos * 2.0 * (1.0 - c) * (1 - _kdamage) +
  //                    _pressure * RankTwoTensor(RankTwoTensor::initIdentity) * c * 2;

  _dstress_dc[_qp] = _pressure * RankTwoTensor(RankTwoTensor::initIdentity) * c * 2;

  RankFourTensor dstress_dstrain;
  // for (unsigned int i = 0; i < _n_phase; ++i)
  //   dstress_dstrain += (*_h_eta[i])[_qp] * (*_dphase_stress_dstrain[i])[_qp];

  dstress_dstrain = _bubble[_qp] * _elasticity_tensor[_qp] * 1.0e-8 +
                    (1 - _bubble[_qp]) * _elasticity_tensor[_qp];

  _Jacobian_mult[_qp] =
      (((1.0 - c) * (1.0 - c) * (1.0 - _kdamage) + _kdamage) * proj_pos + proj_neg) *
      dstress_dstrain;

  _kappa[_qp] = _gc_prop[_qp] * _l[_qp];
  _L[_qp] = 1.0 / (_gc_prop[_qp] * _visco[_qp]);
}
