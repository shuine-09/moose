//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeAnisotropicElasticPFFractureStress.h"
#include "MathUtils.h"

registerMooseObject("TensorMechanicsApp", ComputeAnisotropicElasticPFFractureStress);

template <>
InputParameters
validParams<ComputeAnisotropicElasticPFFractureStress>()
{
  InputParameters params = validParams<ComputeIsotropicLinearElasticPFFractureStress>();
  params.addClassDescription("Computes the stress and free energy derivatives for the phase field "
                             "fracture model, with linear anistropic elasticity");
  return params;
}

ComputeAnisotropicElasticPFFractureStress::ComputeAnisotropicElasticPFFractureStress(
    const InputParameters & parameters)
  : ComputeIsotropicLinearElasticPFFractureStress(parameters)
{
}

void
ComputeAnisotropicElasticPFFractureStress::computeQpStress()
{
  const Real c11 = _elasticity_tensor[_qp](0, 0, 0, 0);
  const Real c12 = _elasticity_tensor[_qp](0, 0, 1, 1);
  const Real c33 = _elasticity_tensor[_qp](2, 2, 2, 2);
  const Real c13 = _elasticity_tensor[_qp](0, 0, 2, 2);

  const Real k0 = ((c11 + c12) * c33 - 2 * c13 * c13) / (c11 + c12 + 2 * c33 - 4 * c13);

  Real c = _c[_qp];

  RankTwoTensor I(RankTwoTensor::initIdentity);
  RankFourTensor I4 = I.outerProduct(I);

  RankFourTensor C0 =
      ((1.0 - c) * (1.0 - c) * (1 - _kdamage) + _kdamage) * _elasticity_tensor[_qp] +
      k0 * I4 * (1 - _kdamage - (1 - _kdamage) * (1 - c) * (1 - c)) *
          MathUtils::heavyside(-_mechanical_strain[_qp].trace());

  _stress[_qp] = C0 * _mechanical_strain[_qp];

  RankFourTensor Ch =
      _elasticity_tensor[_qp] - k0 * I4 * MathUtils::heavyside(-_mechanical_strain[_qp].trace());
  const Real G0_pos =
      (Ch * _mechanical_strain[_qp]).doubleContraction(_mechanical_strain[_qp]) / 2.0;

  // Assign history variable and derivative
  if (G0_pos > _hist_old[_qp])
    _hist[_qp] = G0_pos;
  else
    _hist[_qp] = _hist_old[_qp];

  Real hist_variable = _hist_old[_qp];
  if (_use_current_hist)
    hist_variable = _hist[_qp];

  // Elastic free energy density
  _F[_qp] = hist_variable * ((1.0 - c) * (1.0 - c) * (1 - _kdamage) + _kdamage) +
            _gc[_qp] / (2 * _l[_qp]) * c * c;

  // derivative of elastic free energy density wrt c
  _dFdc[_qp] = -hist_variable * 2.0 * (1.0 - c) * (1 - _kdamage) + _gc[_qp] / _l[_qp] * c;

  // 2nd derivative of elastic free energy density wrt c
  _d2Fdc2[_qp] = hist_variable * 2.0 * (1 - _kdamage) + _gc[_qp] / _l[_qp];

  // 2nd derivative wrt c and strain = 0.0 if we used the previous step's history varible
  if (_use_current_hist)
    _d2Fdcdstrain[_qp] = -(Ch * _mechanical_strain[_qp]) * (1.0 - c) * (1 - _kdamage);

  // Used in StressDivergencePFFracTensors off-diagonal Jacobian
  _dstress_dc[_qp] = (-2 * (1 - c) * (1 - _kdamage) * _elasticity_tensor[_qp] +
                      k0 * I4 * (2 * (1 - c) * (1 - _kdamage)) *
                          MathUtils::heavyside(-_mechanical_strain[_qp].trace())) *
                     _mechanical_strain[_qp];

  _Jacobian_mult[_qp] = C0;
}
