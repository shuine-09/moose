/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "ComputeLinearElasticPFFractureStress.h"

template <>
InputParameters
validParams<ComputeLinearElasticPFFractureStress>()
{
  InputParameters params = validParams<ComputeStressBase>();
  params.addClassDescription("Phase-field fracture model energy contribution to fracture for "
                             "elasticity and undamaged stress under compressive strain");
  params.addRequiredCoupledVar("c", "Order parameter for damage");
  params.addParam<Real>("kdamage", 1e-6, "Stiffness of damaged matrix");
  params.addParam<MaterialPropertyName>(
      "F_name", "E_el", "Name of material property storing the elastic energy");
  return params;
}

ComputeLinearElasticPFFractureStress::ComputeLinearElasticPFFractureStress(
    const InputParameters & parameters)
  : ComputeStressBase(parameters),
    _c(coupledValue("c")),
    _gc_prop(getMaterialProperty<Real>("gc_prop")),

    _l(getMaterialProperty<Real>("l")),

    _kdamage(getParam<Real>("kdamage")),
    _F(declareProperty<Real>(getParam<MaterialPropertyName>("F_name"))),
    _dFdc(declarePropertyDerivative<Real>(getParam<MaterialPropertyName>("F_name"),
                                          getVar("c", 0)->name())),
    _d2Fdc2(declarePropertyDerivative<Real>(
        getParam<MaterialPropertyName>("F_name"), getVar("c", 0)->name(), getVar("c", 0)->name())),
    _d2Fdcdstrain(declareProperty<RankTwoTensor>("d2Fdcdstrain")),
    _dstress_dc(
        declarePropertyDerivative<RankTwoTensor>(_base_name + "stress", getVar("c", 0)->name())),

    _H0_pos(declareProperty<Real>("H0_pos")),
    _H0_pos_old(getMaterialPropertyOld<Real>("H0_pos")),
    _H0_neg(declareProperty<Real>("H0_neg")),
    _H0_neg_old(getMaterialPropertyOld<Real>("H0_neg"))
{
}

void
ComputeLinearElasticPFFractureStress::initQpStatefulProperties()
{
  ComputeStressBase::initQpStatefulProperties();
  _H0_pos[_qp] = 0.0;
}

void
ComputeLinearElasticPFFractureStress::computeQpStress()
{
  Real c = _c[_qp];

  // Compute usual Jacobian mult
  _Jacobian_mult[_qp] = _elasticity_tensor[_qp];

  // Zero out values when c > 1
  Real cfactor = 1.0;
  if (c > 1.0)
    cfactor = 0.0;

  // RankTwoTensor eigvec;
  // std::vector<Real> eigval(LIBMESH_DIM);
  // RankTwoTensor stresscal = _elasticity_tensor[_qp] * _mechanical_strain[_qp];
  // RankTwoTensor I(RankTwoTensor::initIdentity);
  // stresscal.symmetricEigenvaluesEigenvectors(eigval, eigvec);
  // RankTwoTensor Ipos, Ineg, eigval_tensor;
  // eigval_tensor.fillFromInputVector(eigval);
  // // for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  // //   if (eigval[i] >= 0.0)
  // //     Ipos(i, i) = 1.0;
  // //
  // // Ineg = I - Ipos;
  //
  // Ipos.zero();
  // Ineg.zero();
  // for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  // {
  //   if (eigval[i] < 0.0)
  //   {
  //     Ineg(i, i) = 1.0;
  //   }
  // }
  //
  // Ipos = I - Ineg;
  //
  // RankFourTensor QoQ = eigvec.mixedProductIkJl(eigvec);
  // RankFourTensor QToQT = eigvec.transpose().mixedProductIkJl(eigvec.transpose());
  // RankFourTensor IoIp = I.mixedProductIkJl(Ipos);
  // RankFourTensor IoIn = I.mixedProductIkJl(Ineg);
  // // RankFourTensor Jpos = QoQ * IoIp * QToQT * _elasticity_tensor[_qp];
  // RankFourTensor Jneg = QoQ * IoIn * QToQT * _elasticity_tensor[_qp];
  // RankFourTensor Jpos = _elasticity_tensor[_qp] - Jneg;
  //
  // // std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
  // // ((Jpos + Jneg - _elasticity_tensor[_qp]) * RankTwoTensor::initIdentity).print();
  // // std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
  //
  // RankTwoTensor eigpos = eigval_tensor * Ipos;
  // RankTwoTensor eigneg = eigval_tensor * Ineg;
  // // RankTwoTensor stress0pos = eigvec * eigpos * (eigvec.transpose());
  // // RankTwoTensor stress0pos = Jpos * _mechanical_strain[_qp];
  // RankTwoTensor stress0neg = eigvec * eigneg * (eigvec.transpose());
  // RankTwoTensor stress0pos = stresscal - stress0neg;
  // // RankTwoTensor stress0neg = stresscal - stress0pos;
  //
  // Real G0_pos = (stress0pos).doubleContraction(_mechanical_strain[_qp]) / 2.0;
  // Real G0_neg = (stress0neg).doubleContraction(_mechanical_strain[_qp]) / 2.0;

  // std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
  // std::cout << "eigenvalue = " << eigval[0] << ", " << eigval[1] << ", " << eigval[2] <<
  // std::endl; std::cout << "eigenvector = " << std::endl; eigvec.print(); std::cout <<
  // "<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
  //
  // std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
  // std::cout << "stress0pos = " << std::endl;
  // stress0pos.print();
  // std::cout << "proj stress + = " << std::endl;
  // (QoQ * IoIp * QToQT * stresscal).print();
  // std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
  //
  // std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
  // std::cout << "stress0neg = " << std::endl;
  // stress0neg.print();
  // std::cout << "proj stress - = " << std::endl;
  // (QoQ * IoIn * QToQT * stresscal).print();
  // std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
  //
  // std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
  // std::cout << "stress total = " << std::endl;
  // stresscal.print();
  // std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;

  RankTwoTensor eigvec;
  RankFourTensor projpos, projneg;
  RankFourTensor naabb, nabab, nabba;
  RankTwoTensor nab, nba, naa, nbb;
  std::vector<Real> eigval(LIBMESH_DIM);
  std::vector<Real> eigabs(LIBMESH_DIM);
  RankTwoTensor stresscal = _elasticity_tensor[_qp] * _mechanical_strain[_qp];
  stresscal.symmetricEigenvaluesEigenvectors(eigval, eigvec);
  RankTwoTensor Ipos, Ineg, eigval_tensor;
  RankTwoTensor I(RankTwoTensor::initIdentity);
  eigval_tensor.fillFromInputVector(eigval);
  eigabs = eigval;

  Ipos.zero();
  Ineg.zero();
  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  {
    if (eigval[i] < 0.0)
    {
      Ineg(i, i) = 1.0;
      eigabs[i] = -1.0 * eigabs[i];
    }
  }

  Ipos = I - Ineg;

  RankTwoTensor eigpos = eigval_tensor * Ipos;
  RankTwoTensor eigneg = eigval_tensor * Ineg;

  // std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
  // std::cout << "eigpos = " << std::endl;
  // eigpos.print();
  // std::cout << "eigneg = " << std::endl;
  // eigneg.print();
  // std::cout << "total = " << std::endl;
  // eigval_tensor.print();
  // std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;

  projpos.zero();
  projneg.zero();

  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  {
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
    {
      naa.vectorOuterProduct(eigvec.column(i), eigvec.column(i));
      nbb.vectorOuterProduct(eigvec.column(j), eigvec.column(j));
      naabb = naa.outerProduct(nbb);
      projpos += Ipos(i, j) * naabb;
      projneg += Ineg(i, j) * naabb;
    }
  }

  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  {
    for (unsigned int j = 0; j < LIBMESH_DIM && j != i; ++j) // Take care for this part, check it
    {
      // if (MooseUtils::absoluteFuzzyEqual(eigval[i], 0.0) &&
      //     MooseUtils::absoluteFuzzyEqual(eigval[j], 0.0))
      //   continue;
      //
      // Real toldelta = 1e-4; // tolerance value for the eigvalue difference;
      // Real maxeigabs = eigabs[0];
      // for (unsigned int k = 1; k < LIBMESH_DIM; ++k)
      //   {
      //     if (eigabs[k] > maxeigabs)
      //     maxeigabs = eigabs[k];
      //   }
      // if (LIBMESH_DIM ==2 && (abs(eigval[i] - eigval[j])) / maxeigabs < toldelta)
      //   {
      //     eigval_tensor(0,0) *= 1.0 + toldelta;
      //     eigval_tensor(1,1) *= 1.0 - toldelta;
      //   }
      // if (LIBMESH_DIM ==3 && (abs(eigval[i] - eigval[j])) / maxeigabs < toldelta)
      //   {
      //     eigval_tensor(0,0) *= 1.0 + toldelta;
      //     eigval_tensor(1,1) *= 1.0 - toldelta;
      //     eigval_tensor(2,2) = eigval_tensor(2,2)/(1.0 + toldelta)/(1.0 - toldelta);
      //   }
      // if (LIBMESH_DIM == 2 && MooseUtils::absoluteFuzzyEqual(eigval[i], eigval[j]))
      // {
      //   eigval_tensor(0, 0) *= 1.0 + toldelta;
      //   eigval_tensor(1, 1) *= 1.0 - toldelta;
      //   // eigval_tensor(1,1) = eigval_tensor(1,1)/(1.0 + toldelta);
      // }
      // if (LIBMESH_DIM == 3 && MooseUtils::absoluteFuzzyEqual(eigval[i], eigval[j]))
      // {
      //   eigval_tensor(0, 0) *= 1.0 + toldelta;
      //   eigval_tensor(1, 1) *= 1.0 - toldelta;
      //   eigval_tensor(2, 2) = eigval_tensor(2, 2) / (1.0 + toldelta) / (1.0 - toldelta);
      // }
      // if (eigval[i] == eigval[j])
      //   continue;

      nab.vectorOuterProduct(eigvec.column(i), eigvec.column(j));
      nba.vectorOuterProduct(eigvec.column(j), eigvec.column(i));
      // nabab = (nab + nba).outerProduct(nab + nba);
      nabba = nab.outerProduct(nba);
      nabab = (nab).outerProduct(nab + nba);

      if (std::abs(eigval_tensor(i, i) - eigval_tensor(j, j)) > 1.0e-12)
      {
        projpos += (eigpos(i, i) - eigpos(j, j)) / (eigval_tensor(i, i) - eigval_tensor(j, j)) *
                   (nabab) / 2.0;
        projneg += (eigneg(i, i) - eigneg(j, j)) / (eigval_tensor(i, i) - eigval_tensor(j, j)) *
                   (nabab) / 2.0;
      }
      else
      {
        projpos += 1.0 * (nabab) / 2.0;
        projneg += 1.0 * (nabab) / 2.0;
      }
    }
  }

  // RankFourTensor Jpos = projpos * _elasticity_tensor[_qp];
  RankFourTensor Jneg = projneg * _elasticity_tensor[_qp];
  // RankFourTensor Jneg = _elasticity_tensor[_qp] - Jpos;
  RankFourTensor Jpos = _elasticity_tensor[_qp] - Jneg;

  // std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
  // ((projpos + projneg - RankFourTensor::initIdentity) * stresscal).print();
  // std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;

  // RankTwoTensor eigpos = eigval_tensor * Ipos;
  // RankTwoTensor eigneg = eigval_tensor * Ineg;
  // RankTwoTensor stresspos = eigvec * eigpos * (eigvec.transpose()); // principle positive
  // stress RankTwoTensor stressneg = eigvec * eigneg * (eigvec.transpose()); // principle
  // negative stress RankTwoTensor stress0pos = eigvec * eigpos * (eigvec.transpose()); //Postive
  // part of strain == Q*eigpos*QT RankTwoTensor stress0neg = eigvec * eigneg *
  // (eigvec.transpose()); //Negative part of strain == Q*eigpos*QT RankTwoTensor stress0pos =
  // (Jpos * _mechanical_strain[_qp] + stresspos) / 2.0; RankTwoTensor stress0neg = (Jneg *
  // _mechanical_strain[_qp] + stressneg) / 2.0;

  // RankTwoTensor stress0pos = eigvec * eigpos * (eigvec.transpose());
  // RankTwoTensor stress0neg = stresscal - stress0pos;
  RankTwoTensor stress0neg = eigvec * eigneg * (eigvec.transpose());
  RankTwoTensor stress0pos = stresscal - stress0neg;

  // std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
  // std::cout << "stress0pos = " << std::endl;
  // stress0pos.print();
  // std::cout << "proj stress + = " << std::endl;
  // (projpos * stresscal).print();
  // std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
  //
  // std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
  // std::cout << "stress0neg = " << std::endl;
  // stress0neg.print();
  // std::cout << "proj stress - = " << std::endl;
  // (projneg * stresscal).print();
  // std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
  //
  // std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
  // std::cout << "stress total = " << std::endl;
  // stresscal.print();
  // std::cout << "proj stress total = " << std::endl;
  // (projneg * stresscal + projpos * stresscal).print();
  // std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;

  // std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
  // std::cout << "stress plus = " << std::endl;
  // (Jpos * _mechanical_strain[_qp]).print();
  // std::cout << "stress minus = " << std::endl;
  // (Jneg * _mechanical_strain[_qp]).print();
  // std::cout << "stress total = " << std::endl;
  // stresscal.print();
  // std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;

  Real G0_pos = (stress0pos).doubleContraction(_mechanical_strain[_qp]) / 2.0;
  Real G0_neg = (stress0neg).doubleContraction(_mechanical_strain[_qp]) / 2.0;

  // if (c < 1.0)
  _Jacobian_mult[_qp] = ((1.0 - c) * (1.0 - c) * (1.0 - _kdamage) + _kdamage) * Jpos + Jneg;
  // else
  //   _Jacobian_mult[_qp] = _elasticity_tensor[_qp];

  // _Jacobian_mult[_qp] =
  //     (cfactor * (1.0 - c) * (1.0 - c) * (1.0 - _kdamage) + _kdamage) * (Jpos + Jneg);

  // _H0_neg[_qp] = G0_neg;
  //
  // if (_H0_pos_old[_qp] < _H0_neg_old[_qp])
  //   c = 0.0;

  // _Jacobian_mult[_qp] =
  //     ((1.0 - c) * (1.0 - c) * (1.0 - _kdamage) + _kdamage) * _elasticity_tensor[_qp];

  if (G0_pos > _H0_pos_old[_qp])
    _H0_pos[_qp] = G0_pos;
  else
    _H0_pos[_qp] = _H0_pos_old[_qp];

  // _H0_pos[_qp] = stresscal.doubleContraction(_mechanical_strain[_qp]) / 2.0;

  // Damage associated with positive component of stress
  _stress[_qp] = stress0pos * (cfactor * (1.0 - c) * (1.0 - c) + _kdamage) + stress0neg;
  // _stress[_qp] = (cfactor * (1.0 - c) * (1.0 - c) * (1.0 - _kdamage) + _kdamage) * stresscal;

  // Elastic free energy density
  _F[_qp] = _H0_pos_old[_qp] * (cfactor * (1.0 - c) * (1.0 - c) + _kdamage) + G0_neg +
            _gc_prop[_qp] * c * c / 2.0 / _l[_qp];

  // derivative of elastic free energy density wrt c
  _dFdc[_qp] = -_H0_pos_old[_qp] * 2.0 * cfactor * (1.0 - c) + _gc_prop[_qp] * c / _l[_qp];

  // 2nd derivative of elastic free energy density wrt c
  _d2Fdc2[_qp] = _H0_pos_old[_qp] * 2.0 * cfactor + _gc_prop[_qp] / _l[_qp];

  // 2nd derivative wrt c and strain
  if (G0_pos > _H0_pos_old[_qp])
    //_d2Fdcdstrain[_qp] = -2.0 * stress0pos * cfactor * (1.0 - c) * (1.0 - _kdamage);
    _d2Fdcdstrain[_qp] = stress0pos * 0.5;
  else
    _d2Fdcdstrain[_qp] = stress0pos * 0.0;

  // Used in StressDivergencePFFracTensors off-diagonal Jacobian
  _dstress_dc[_qp] = -stress0pos * 2.0 * cfactor * (1.0 - c);
  // _dstress_dc[_qp] = -stresscal * 2.0 * cfactor * (1.0 - c) * (1.0 - _kdamage);
}
