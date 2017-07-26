/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "AnisoStressElasticPFDamage.h"
#include "libmesh/utility.h"

template <>
InputParameters
validParams<AnisoStressElasticPFDamage>()
{
  InputParameters params = validParams<ComputeStressBase>();
  params.addClassDescription("Phase-field fracture model energy contribution to damage "
                             "growth-isotropic elasticity and undamaged stress under compressive "
                             "strain");
  params.addRequiredCoupledVar("c", "Order parameter for damage");
  params.addParam<Real>("kdamage", 1e-6, "Stiffness of damaged matrix");

  return params;
}

AnisoStressElasticPFDamage::AnisoStressElasticPFDamage(const InputParameters & parameters)
  : ComputeStressBase(parameters),
    _c(coupledValue("c")),
    _kdamage(getParam<Real>("kdamage")),
    _G0_pos(declareProperty<Real>("G0_pos")),
    _G0_pos_old(declarePropertyOld<Real>("G0_pos")), // added G0_old variable for old time step
    _dstress_dc(
        declarePropertyDerivative<RankTwoTensor>(_base_name + "stress", getVar("c", 0)->name())),
    _dG0_pos_dstrain(declareProperty<RankTwoTensor>("dG0_pos_dstrain")),
    _etens(LIBMESH_DIM),
    _epos(LIBMESH_DIM),
    _eigval(LIBMESH_DIM)
{
}

void
AnisoStressElasticPFDamage::initQpStatefulProperties()
{
  ComputeStressBase::initQpStatefulProperties();
  // if (_t > 0)
  //_G0_pos_old[_qp] = _G0_pos[_qp];
  // else{
  _G0_pos[_qp] = 0.0;
  _G0_pos_old[_qp] = 0.0;
  //}
  // _dG0_pos_dstrain[_qp] = _stress[_qp];
  // _dstress_dc[_qp] = -_stress[_qp] * (2.0 * (1.0 - _c[_qp]));
}

void
AnisoStressElasticPFDamage::computeQpStress()
{
  updateVar();
  updateJacobian();
}

void
AnisoStressElasticPFDamage::updateVar()
{

  // Real c = _c[_qp];
  // // Real cfactor = 1.0;
  // // if (c > 1.0)
  // //   cfactor = 0.0;
  //
  // // _Jacobian_mult[_qp] = _elasticity_tensor[_qp];
  //
  // RankTwoTensor eigvec;
  // // RankFourTensor projpos, projneg, projpos1, projpos2, projneg1, projneg2;
  // // RankFourTensor naabb, nabab, nabba;
  // // RankTwoTensor nab, nba, naa, nbb;
  // std::vector<Real> eigval(LIBMESH_DIM);
  // std::vector<Real> eigabs(LIBMESH_DIM);
  // RankTwoTensor stresscal = _elasticity_tensor[_qp] * _mechanical_strain[_qp];
  // stresscal.symmetricEigenvaluesEigenvectors(eigval, eigvec);
  // RankTwoTensor Ipos, Ineg, eigval_tensor;
  // RankTwoTensor I(RankTwoTensor::initIdentity);
  // eigval_tensor.fillFromInputVector(eigval);
  // eigabs = eigval;
  // for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  //   {
  //     if (eigval[i] < 0.0)
  //     {
  //       Ineg(i, i) = 1.0;
  //     }
  //   }
  //
  // Ipos = I - Ineg;
  //
  //
  //
  // RankFourTensor QoQ = eigvec.mixedProductIkJl(eigvec);
  // RankFourTensor QToQT = eigvec.transpose().mixedProductIkJl(eigvec.transpose());
  // RankFourTensor IpoIp = I.mixedProductIkJl(Ipos);
  // RankFourTensor InoIn = I.mixedProductIkJl(Ineg);
  //   // RankFourTensor IpIp = Ipos.mixedProductIkJl(I);
  //   // RankFourTensor InIn = Ipos.mixedProductIkJl(I);
  //   // RankFourTensor Jpos = QoQ * (IpoIp + IpIp) * QToQT * _elasticity_tensor[_qp] / 2.0;
  //   // RankFourTensor Jneg = QoQ * (InoIn + InIn) * QToQT * _elasticity_tensor[_qp] / 2.0;
  //   // RankFourTensor Jpos = QoQ * IpIp * QToQT * _elasticity_tensor[_qp];
  //   // RankFourTensor Jneg = QoQ * InIn * QToQT * _elasticity_tensor[_qp];
  // RankFourTensor Jpos = QoQ * IpoIp * QToQT * _elasticity_tensor[_qp];
  // RankFourTensor Jneg = QoQ * InoIn * QToQT * _elasticity_tensor[_qp];

  Real c = _c[_qp];
  // Real cfactor = 1.0;

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
  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  {
    if (eigval[i] < 0.0)
    {
      Ineg(i, i) = 1.0;
      eigabs[i] = -1.0 * eigabs[i];
    }
  }

  Ipos = I - Ineg;

  // RankFourTensor QoQ = eigvec.mixedProductIkJl(eigvec);
  // RankFourTensor QToQT = eigvec.transpose().mixedProductIkJl(eigvec.transpose());
  // RankFourTensor IpoIp = I.mixedProductIkJl(Ipos);
  // RankFourTensor InoIn = I.mixedProductIkJl(Ineg);
  // // RankFourTensor IpIp = Ipos.mixedProductIkJl(I);
  // // RankFourTensor InIn = Ipos.mixedProductIkJl(I);
  // // RankFourTensor Jpos = QoQ * (IpoIp + IpIp) * QToQT * _elasticity_tensor[_qp] / 2.0;
  // // RankFourTensor Jneg = QoQ * (InoIn + InIn) * QToQT * _elasticity_tensor[_qp] / 2.0;
  // // RankFourTensor Jpos = QoQ * IpIp * QToQT * _elasticity_tensor[_qp];
  // // RankFourTensor Jneg = QoQ * InIn * QToQT * _elasticity_tensor[_qp];
  // RankFourTensor Jpos = QoQ * IpoIp * QToQT * _elasticity_tensor[_qp];
  // RankFourTensor Jneg = QoQ * InoIn * QToQT * _elasticity_tensor[_qp];

  RankTwoTensor eigpos = eigval_tensor * Ipos;
  RankTwoTensor eigneg = eigval_tensor * Ineg;

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
      if (MooseUtils::absoluteFuzzyEqual(eigval[i], 0.0) &&
          MooseUtils::absoluteFuzzyEqual(eigval[j], 0.0))
        continue;

      Real toldelta = 1e-4; // tolerance value for the eigvalue difference;
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
      if (LIBMESH_DIM == 2 && MooseUtils::absoluteFuzzyEqual(eigval[i], eigval[j]))
      {
        eigval_tensor(0, 0) *= 1.0 + toldelta;
        eigval_tensor(1, 1) *= 1.0 - toldelta;
        // eigval_tensor(1,1) = eigval_tensor(1,1)/(1.0 + toldelta);
      }
      if (LIBMESH_DIM == 3 && MooseUtils::absoluteFuzzyEqual(eigval[i], eigval[j]))
      {
        eigval_tensor(0, 0) *= 1.0 + toldelta;
        eigval_tensor(1, 1) *= 1.0 - toldelta;
        eigval_tensor(2, 2) = eigval_tensor(2, 2) / (1.0 + toldelta) / (1.0 - toldelta);
      }
      // if (eigval[i] == eigval[j])
      //   continue;

      nab.vectorOuterProduct(eigvec.column(i), eigvec.column(j));
      nba.vectorOuterProduct(eigvec.column(j), eigvec.column(i));
      // nabab = (nab + nba).outerProduct(nab + nba);
      nabba = nab.outerProduct(nba);

      projpos += (eigpos(i, i) - eigpos(j, j)) / (eigval_tensor(i, i) - eigval_tensor(j, j)) *
                 (nabab + nabba) / 2.0;
      projneg += (eigneg(i, i) - eigneg(j, j)) / (eigval_tensor(i, i) - eigval_tensor(j, j)) *
                 (nabab + nabba) / 2.0;
    }
  }
  RankFourTensor Jpos = projpos * _elasticity_tensor[_qp];
  RankFourTensor Jneg = projneg * _elasticity_tensor[_qp];

  // RankTwoTensor eigpos = eigval_tensor * Ipos;
  // RankTwoTensor eigneg = eigval_tensor * Ineg;
  RankTwoTensor stresspos = eigvec * eigpos * (eigvec.transpose()); // principle positive stress
  RankTwoTensor stressneg = eigvec * eigneg * (eigvec.transpose()); // principle negative stress
  // RankTwoTensor stress0pos = eigvec * eigpos * (eigvec.transpose()); //Postive part of strain ==
  // Q*eigpos*QT RankTwoTensor stress0neg = eigvec * eigneg * (eigvec.transpose()); //Negative part
  // of strain == Q*eigpos*QT
  RankTwoTensor stress0pos = (Jpos * _mechanical_strain[_qp] + stresspos) / 2.0;
  RankTwoTensor stress0neg = (Jneg * _mechanical_strain[_qp] + stressneg) / 2.0;

  Real G0_Pos_trial = (stresspos).doubleContraction(_mechanical_strain[_qp]) / 2.0;

  if (G0_Pos_trial > _G0_pos_old[_qp])
  {
    _G0_pos[_qp] = G0_Pos_trial;
    _dG0_pos_dstrain[_qp] = stress0pos;
  }
  else
  {
    _G0_pos[_qp] = _G0_pos_old[_qp];
    _dG0_pos_dstrain[_qp] = stress0pos * 0.0; // derivative of positive strain energy wrt strain
  }

  // Used in PFFracBulkRate Jacobian
  // _dG0_pos_dstrain[_qp] = stress0pos;

  // Damage associated with positive component of stress
  _stress[_qp] = stress0pos * ((1.0 - c) * (1.0 - c) * (1.0 - _kdamage) + _kdamage) + stress0neg;

  // Used in StressDivergencePFFracTensors Jacobian
  if (c < 1.0)
    _dstress_dc[_qp] = -stress0pos * 2.0 * (1.0 - c) * (1.0 - _kdamage);
  else
    _dstress_dc[_qp].zero();
  // _dstress_dc[_qp] = -cfactor * _dG0_pos_dstrain[_qp] * 2.0 * (1.0 - c);
}

void
AnisoStressElasticPFDamage::updateJacobian()
{
  Real c = _c[_qp];
  // Real cfactor = 1.0;

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
  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  {
    if (eigval[i] < 0.0)
    {
      Ineg(i, i) = 1.0;
      eigabs[i] = -1.0 * eigabs[i];
    }
  }

  Ipos = I - Ineg;

  RankFourTensor QoQ = eigvec.mixedProductIkJl(eigvec);
  RankFourTensor QToQT = eigvec.transpose().mixedProductIkJl(eigvec.transpose());
  RankFourTensor IpoIp = I.mixedProductIkJl(Ipos);
  RankFourTensor InoIn = I.mixedProductIkJl(Ineg);
  // RankFourTensor IpIp = Ipos.mixedProductIkJl(I);
  // RankFourTensor InIn = Ipos.mixedProductIkJl(I);
  // RankFourTensor Jpos = QoQ * (IpoIp + IpIp) * QToQT * _elasticity_tensor[_qp] / 2.0;
  // RankFourTensor Jneg = QoQ * (InoIn + InIn) * QToQT * _elasticity_tensor[_qp] / 2.0;
  // RankFourTensor Jpos = QoQ * IpIp * QToQT * _elasticity_tensor[_qp];
  // RankFourTensor Jneg = QoQ * InIn * QToQT * _elasticity_tensor[_qp];
  RankFourTensor Jpos = QoQ * IpoIp * QToQT * _elasticity_tensor[_qp];
  RankFourTensor Jneg = QoQ * InoIn * QToQT * _elasticity_tensor[_qp];

  // RankTwoTensor eigpos = eigval_tensor * Ipos;
  // RankTwoTensor eigneg = eigval_tensor * Ineg;
  //
  // for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  //   {
  //      for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
  //        {
  //          naa.vectorOuterProduct(eigvec.column(i), eigvec.column(i));
  //          nbb.vectorOuterProduct(eigvec.column(j), eigvec.column(j));
  //          naabb = naa.outerProduct(nbb);
  //          projpos += Ipos(i,j) * naabb;
  //          projneg += Ineg(i,j) * naabb;
  //         }
  //   }
  //
  //
  //
  //
  // for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  //   {
  //     for (unsigned int j = 0; j < LIBMESH_DIM && j!= i; ++j) //Take care for this part, check it
  //       {
  //         if (MooseUtils::absoluteFuzzyEqual(eigval[i], 0.0) &&
  //         MooseUtils::absoluteFuzzyEqual(eigval[j], 0.0)) continue;
  //
  //         Real toldelta = 1e-4; //tolerance value for the eigvalue difference;
  //         // Real maxeigabs = eigabs[0];
  //         // for (unsigned int k = 1; k < LIBMESH_DIM; ++k)
  //         //   {
  //         //     if (eigabs[k] > maxeigabs)
  //         //     maxeigabs = eigabs[k];
  //         //   }
  //         // if (LIBMESH_DIM ==2 && (abs(eigval[i] - eigval[j])) / maxeigabs < toldelta)
  //         //   {
  //         //     eigval_tensor(0,0) *= 1.0 + toldelta;
  //         //     eigval_tensor(1,1) *= 1.0 - toldelta;
  //         //   }
  //         // if (LIBMESH_DIM ==3 && (abs(eigval[i] - eigval[j])) / maxeigabs < toldelta)
  //         //   {
  //         //     eigval_tensor(0,0) *= 1.0 + toldelta;
  //         //     eigval_tensor(1,1) *= 1.0 - toldelta;
  //         //     eigval_tensor(2,2) = eigval_tensor(2,2)/(1.0 + toldelta)/(1.0 - toldelta);
  //         //   }
  //         if (LIBMESH_DIM ==2 && MooseUtils::absoluteFuzzyEqual(eigval[i], eigval[j]))
  //           {
  //             eigval_tensor(0,0) *= 1.0 + toldelta;
  //             eigval_tensor(1,1) *= 1.0 - toldelta;
  //             // eigval_tensor(1,1) = eigval_tensor(1,1)/(1.0 + toldelta);
  //           }
  //         if (LIBMESH_DIM ==3 && MooseUtils::absoluteFuzzyEqual(eigval[i], eigval[j]))
  //           {
  //             eigval_tensor(0,0) *= 1.0 + toldelta;
  //             eigval_tensor(1,1) *= 1.0 - toldelta;
  //             eigval_tensor(2,2) = eigval_tensor(2,2)/(1.0 + toldelta)/(1.0 - toldelta);
  //           }
  //               // if (eigval[i] == eigval[j])
  //               //   continue;
  //
  //         nab.vectorOuterProduct(eigvec.column(i), eigvec.column(j));
  //         nba.vectorOuterProduct(eigvec.column(j), eigvec.column(i));
  //                   // nabab = (nab + nba).outerProduct(nab + nba);
  //         nabba = nab.outerProduct(nba);
  //
  //         projpos += (eigpos(i,i) - eigpos(j,j))/(eigval_tensor(i,i) - eigval_tensor(j,j)) *
  //         (nabab + nabba) / 2.0; projneg += (eigneg(i,i) - eigneg(j,j))/(eigval_tensor(i,i) -
  //         eigval_tensor(j,j)) * (nabab + nabba) / 2.0;
  //       }
  //   }
  // RankFourTensor Jpos = projpos * _elasticity_tensor[_qp];
  // RankFourTensor Jneg = projneg * _elasticity_tensor[_qp];

  if (c < 1.0)
    _Jacobian_mult[_qp] = ((1.0 - c) * (1.0 - c) * (1.0 - _kdamage) + _kdamage) * Jpos + Jneg;
  else
    _Jacobian_mult[_qp] = _elasticity_tensor[_qp];
}
