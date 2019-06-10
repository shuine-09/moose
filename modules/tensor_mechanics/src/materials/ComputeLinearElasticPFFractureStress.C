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
  InputParameters params = validParams<ComputeStressBase>();
  params.addClassDescription("Phase-field fracture model energy contribution to fracture for "
                             "elasticity and undamaged stress under compressive strain");
  params.addRequiredCoupledVar("c", "Order parameter for damage");
  params.addParam<Real>("kdamage", 1e-6, "Stiffness of damaged matrix");
  params.addParam<bool>(
      "use_current_history_variable", true, "Use the current value of the history variable.");
  params.addParam<bool>(
      "decomp", true, "Use the current value of the history variable.");
  params.addParam<bool>(
      "linear", false, "Use linear fracture energy.");
  params.addParam<bool>(
      "use_vi", false, "Use vi solver.");
  params.addParam<MaterialPropertyName>(
      "F_name", "E_el", "Name of material property storing the elastic energy");
  params.addParam<MaterialPropertyName>(
      "kappa_name",
      "kappa_op",
      "Name of material property being created to store the interfacial parameter kappa");
  params.addParam<MaterialPropertyName>(
      "mobility_name", "L", "Name of material property being created to store the mobility L");
  return params;
}

ComputeLinearElasticPFFractureStress::ComputeLinearElasticPFFractureStress(
    const InputParameters & parameters)
  : ComputeStressBase(parameters),
    _use_current_hist(getParam<bool>("use_current_history_variable")),
    _decomp(getParam<bool>("decomp")),
     _linear(getParam<bool>("linear")),
    _use_vi(getParam<bool>("use_vi")),
    _c(coupledValue("c")),
    _gc_prop(getMaterialProperty<Real>("gc_prop")),
    _l(getMaterialProperty<Real>("l")),
    _visco(getMaterialProperty<Real>("visco")),
    _kdamage(getParam<Real>("kdamage")),
    _F(declareProperty<Real>(getParam<MaterialPropertyName>("F_name"))),
    _dFdc(declarePropertyDerivative<Real>(getParam<MaterialPropertyName>("F_name"),
                                          getVar("c", 0)->name())),
    _d2Fdc2(declarePropertyDerivative<Real>(
        getParam<MaterialPropertyName>("F_name"), getVar("c", 0)->name(), getVar("c", 0)->name())),
    _d2Fdcdstrain(declareProperty<RankTwoTensor>("d2Fdcdstrain")),
    _dstress_dc(declarePropertyDerivative<RankTwoTensor>("stress", getVar("c", 0)->name())),
    _hist(declareProperty<Real>("hist")),
    _hist_old(getMaterialPropertyOld<Real>("hist")),
    _kappa(declareProperty<Real>(getParam<MaterialPropertyName>("kappa_name"))),
    _L(declareProperty<Real>(getParam<MaterialPropertyName>("mobility_name"))),
    _proj(declareProperty<RankFourTensor>("proj")),
    _proj_old(getMaterialPropertyOld<RankFourTensor>("proj")),
    _proj_older(getMaterialPropertyOlder<RankFourTensor>("proj"))
{
}

void
ComputeLinearElasticPFFractureStress::computeQpStress()
{
  const Real c = _c[_qp];

  // Zero out values when c > 1
  Real cfactor = 1.0;
  if (c > 1.0)
    cfactor = 0.0;

  // Compute Uncracked stress
  RankTwoTensor uncracked_stress = _elasticity_tensor[_qp] * _mechanical_strain[_qp];

  // Create the positive and negative projection tensors
  RankFourTensor I4sym(RankFourTensor::initIdentitySymmetricFour);
  std::vector<Real> eigval;
  RankTwoTensor eigvec;
  RankFourTensor Ppos = uncracked_stress.positiveProjectionEigenDecomposition(eigval, eigvec);
  RankFourTensor Pneg = I4sym - Ppos;

  _proj[_qp] = Ppos;

  //if (c > 0.9)
  //  Ppos = _proj_old[_qp];

  Pneg = I4sym - Ppos;

  // Project the positive and negative stresses
  RankTwoTensor stress0pos = Ppos * uncracked_stress;
  RankTwoTensor stress0neg = Pneg * uncracked_stress;

  // Compute the positive and negative elastic energies
  Real G0_pos = (stress0pos).doubleContraction(_mechanical_strain[_qp]) / 2.0;
  Real G0_neg = (stress0neg).doubleContraction(_mechanical_strain[_qp]) / 2.0;

  if (_use_vi)
    _hist[_qp] = G0_pos;
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

  // Compute degradation function and derivatives
  Real h = cfactor * (1.0 - c) * (1.0 - c) * (1.0 - _kdamage) + _kdamage;
  Real dhdc = -2.0 * cfactor * (1.0 - c) * (1.0 - _kdamage);
  Real d2hdc2 = 2.0 * cfactor * (1.0 - _kdamage);

  // Compute stress and its derivatives
  if (_decomp)
  {
    _stress[_qp] = stress0pos * h + stress0neg; // equivalent to (Ppos * h + Pneg) * uncracked_stress;
    _dstress_dc[_qp] = stress0pos * dhdc;
    _Jacobian_mult[_qp] = (Ppos * h + Pneg) * _elasticity_tensor[_qp];

  }
  else
  {
    _stress[_qp] = (stress0pos + stress0neg) * h; // equivalent to (Ppos * h + Pneg) * uncracked_stress;
    _dstress_dc[_qp] = _stress[_qp] * dhdc;
    _Jacobian_mult[_qp] = h * _elasticity_tensor[_qp];
  }

  // Compute energy and its derivatives
  if (_linear)
  {
    _F[_qp] = hist_variable * h - G0_neg + 3 * _gc_prop[_qp] * c / (8 * _l[_qp]);
    _dFdc[_qp] = hist_variable * dhdc + 3 * _gc_prop[_qp] / (8*_l[_qp]);
    _d2Fdc2[_qp] = hist_variable * d2hdc2;
  }
  else
  {  _F[_qp] = hist_variable * h - G0_neg + _gc_prop[_qp] * c * c / (2 * _l[_qp]);
    _dFdc[_qp] = hist_variable * dhdc + _gc_prop[_qp] * c / _l[_qp];
    _d2Fdc2[_qp] = hist_variable * d2hdc2 + _gc_prop[_qp] / _l[_qp];
  }

  // 2nd derivative wrt c and strain = 0.0 if we used the previous step's history varible
  if (_use_current_hist)
    _d2Fdcdstrain[_qp] = _dstress_dc[_qp];

  // Assign L and kappa
  _kappa[_qp] = _gc_prop[_qp] * _l[_qp];
  _L[_qp] = 1.0 / (_gc_prop[_qp] * _visco[_qp]);
}
