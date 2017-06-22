/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "AnisoElasticPFDamage.h"
#include "libmesh/utility.h"

template <>
InputParameters
validParams<AnisoElasticPFDamage>()
{
  InputParameters params = validParams<ComputeStressBase>();
  params.addClassDescription("Phase-field fracture model energy contribution to damage "
                             "growth-isotropic elasticity and undamaged stress under compressive "
                             "strain");
  params.addRequiredCoupledVar("c", "Order parameter for damage");
  params.addParam<Real>("kdamage", 1e-6, "Stiffness of damaged matrix");

  return params;
}

AnisoElasticPFDamage::AnisoElasticPFDamage(const InputParameters & parameters)
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
AnisoElasticPFDamage::initQpStatefulProperties()
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
AnisoElasticPFDamage::computeQpStress()
{
  updateVar();
  updateJacobian();
}

void
AnisoElasticPFDamage::updateVar()
{

  Real c = _c[_qp];
  Real cfactor = 1.0;
  // if (c > 1.0)
  //   cfactor = 0.0;

  // _Jacobian_mult[_qp] = _elasticity_tensor[_qp];

  RankTwoTensor eigvec;
  RankFourTensor projpos, projneg;
  RankFourTensor naabb, nabab, nabba;
  RankTwoTensor nab, nba, naa, nbb;
  std::vector<Real> eigval(LIBMESH_DIM);
  std::vector<Real> eigabs(LIBMESH_DIM);
  _mechanical_strain[_qp].symmetricEigenvaluesEigenvectors(eigval, eigvec);
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
  // Creation of some critical fourth order tensors
  // RankFourTensor QoQ = eigvec.mixedProductIkJl(eigvec);
  // RankFourTensor QToQT = (eigvec.transpose()).mixedProductIkJl(eigvec.transpose());
  // // RankFourTensor QToQ = (eigvec.transpose()).mixedProductIkJl(eigvec);
  // RankFourTensor IpoIp = Ipos.mixedProductIkJl(I);
  // RankFourTensor InoIn = Ineg.mixedProductIkJl(I);
  // RankFourTensor proj_pos = QoQ * IpoIp * QToQT;  //positive projection tensor -- Fourthorder
  // //Question about IpoIp, which may should be cross product RankFourTensor proj_neg = QoQ * InoIn
  // * QToQT;  //Negative projection tensor -- Fourthorder //InoIn may should be
  // crossproduct/Inneproduct

  RankFourTensor IooI = I.mixedProductIkJl(I);
  // RankFourTensor IoxI = I.outerProduct(I);

  RankTwoTensor eigpos = eigval_tensor * Ipos;
  RankTwoTensor eigneg = eigval_tensor * Ineg;
  // RankTwoTensor strainpos = QoQ * _mechanical_strain[_qp];
  RankTwoTensor strainpos =
      eigvec * eigpos * (eigvec.transpose()); // Postive part of strain == Q*eigpos*QT
  RankTwoTensor strainneg =
      eigvec * eigneg * (eigvec.transpose()); // Negative part of strain == Q*eigpos*QT

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  // Calculate the fourth-order tensor projection tensor
  ////////////////////////////////////////////////////////////////
  // for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  //     {
  //       for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
  //         {
  //           nab.fourvectorOuterProduct(eigvec.column(i),eigvec.column(i),eigvec.column(j),eigvec.column(j));
  //           projpos += Ipos(i,j) * nab;
  //           projneg += Ineg(i,j) * nab;
  //         }
  //       for (unsigned int j = 0 && j !=i; j < LIBMESH_DIM; ++j)
  //         {
  //           nba.fourvectorOuterProduct(eigvec.column(i),eigvec.column(j),eigvec.column(i),eigvec.column(j));
  //           nbba.fourvectorOuterProduct(eigvec.column(i),eigvec.column(j),eigvec.column(j),eigvec.column(i));
  //           projpos += (Ipos(i,i) - Ipos(j,j))/(eigval_tensor(i,i) - eigval_tensor(j,j)) * (nba +
  //           nbba) / 2.0; projneg += (Ineg(i,i) - Ineg(j,j))/(eigval_tensor(i,i) -
  //           eigval_tensor(j,j)) * (nba + nbba) / 2.0;
  //         }
  //     }
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  // for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  //     {
  //       for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
  //         {
  //           for (unsigned int p = 0; p < LIBMESH_DIM; ++p)
  //             {
  //               for (unsigned int q = 0; q < LIBMESH_DIM; ++q)
  //                 {
  //                   for (unsigned int m = 0; m < LIBMESH_DIM; ++m)
  //                     {
  //                       for (unsigned int n = 0; n < LIBMESH_DIM; ++n)
  //                         {
  //                           naabb(p, q, m, n) += eigvec.column(i)(p) * eigvec.column(i)(q) *
  //                           eigvec.column(j)(m) * eigvec.column(j)(n);
  //                         }
  //                     }
  //                  }
  //               }
  //           projpos += Ipos(i,j) * naabb;
  //           projneg += Ineg(i,j) * naabb;
  //         }
  //     }
  //
  //   for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  //     {
  //       for (unsigned int j = 0; j < LIBMESH_DIM  && j !=i; ++j) //Take care for this part, check
  //       it
  //
  //         {
  //           nab.vectorOuterProduct(eigvec.column(i), eigvec.column(j));
  //           nba.vectorOuterProduct(eigvec.column(j), eigvec.column(i));
  //           nabab = nab.outerProduct(nab);
  //           nabba = nab.outerProduct(nba);
  //
  //           // nba.fourvectorOuterProduct(eigvec.column(i),eigvec.column(j),eigvec.column(i),eigvec.column(j));
  //           // nbba.fourvectorOuterProduct(eigvec.column(i),eigvec.column(j),eigvec.column(j),eigvec.column(i));
  //           // projpos += (Ipos(i,i) * eigval_tensor(i,i) - Ipos(j,j) * eigval_tensor(j,j))/(eigval_tensor(i,i) - eigval_tensor(j,j)) * (nabab + nabba) / 2.0;
  //           // projneg += (Ineg(i,i) * eigval_tensor(i,i) - Ineg(j,j) * eigval_tensor(j,j))/(eigval_tensor(i,i) - eigval_tensor(j,j)) * (nabab + nabba) / 2.0;
  //           projpos += (eigpos(i,i) - eigpos(j,j))/(eigval_tensor(i,i) - eigval_tensor(j,j)) *
  //           (nabab + nabba) / 2.0; projneg += (eigneg(i,i) - eigneg(j,j))/(eigval_tensor(i,i) -
  //           eigval_tensor(j,j)) * (nabab + nabba) / 2.0;
  //
  //         }
  //     }

  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////

  // for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  //     {
  //       for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
  //         {
  //           naa.vectorOuterProduct(eigvec.column(i), eigvec.column(i));
  //           nbb.vectorOuterProduct(eigvec.column(j), eigvec.column(j));
  //           naabb = naa.outerProduct(nbb);
  //           projpos += Ipos(i,j) * naabb;
  //           projneg += Ineg(i,j) * naabb;
  //         }
  //     }
  //
  // for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  //   {
  //      for (unsigned int j = 0; j < LIBMESH_DIM; ++j) //Take care for this part, check it
  //         {
  //           if (eigval[i] == eigval[j])
  //             continue;
  //
  //           nab.vectorOuterProduct(eigvec.column(i), eigvec.column(j));
  //           nba.vectorOuterProduct(eigvec.column(j), eigvec.column(i));
  //           nabba = nab.outerProduct(nba);
  //           projpos += (eigpos(i,i) - eigpos(j,j))/(eigval_tensor(i,i) - eigval_tensor(j,j)) *
  //           (nabab + nabba) / 2.0; projneg += (eigneg(i,i) - eigneg(j,j))/(eigval_tensor(i,i) -
  //           eigval_tensor(j,j)) * (nabab + nabba) / 2.0;
  //       }
  //   }

  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////

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

      Real toldelta = 1e-6; // tolerance value for the eigvalue difference;
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
  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////

  // Define positive and negative part of stress
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  // RankTwoTensor stress0pos = ((_elasticity_tensor[_qp] * strainpos).initialContraction(proj_pos)
  // )/2.0; RankTwoTensor stress0neg = (_elasticity_tensor[_qp] * proj_pos * strainneg +
  // (_elasticity_tensor[_qp] * strainpos).initialContraction(proj_neg) + _elasticity_tensor[_qp] *
  // proj_neg * strainpos + (_elasticity_tensor[_qp] * strainneg).initialContraction(proj_pos) +
  // (_elasticity_tensor[_qp] * strainneg).initialContraction(proj_neg))/2.0;

  ////////////////////////////////////////////////////////////////////////////////
  // RankTwoTensor stress0pos = ((_elasticity_tensor[_qp] *
  // projpos).chainruleContraction(_mechanical_strain[_qp]) + (_elasticity_tensor[_qp] *
  // strainpos).inoutContraction(IooI))/2.0; RankTwoTensor stress0neg = ((_elasticity_tensor[_qp] *
  // projneg).chainruleContraction(_mechanical_strain[_qp]) + (_elasticity_tensor[_qp] *
  // strainneg).inoutContraction(IooI))/2.0;
  RankTwoTensor stress0pos =
      (_mechanical_strain[_qp].initialContraction(_elasticity_tensor[_qp] * projpos) +
       (_elasticity_tensor[_qp] * strainpos).initialContraction(IooI)) /
      2.0;
  RankTwoTensor stress0neg =
      (_mechanical_strain[_qp].initialContraction(_elasticity_tensor[_qp] * projneg) +
       (_elasticity_tensor[_qp] * strainneg).initialContraction(IooI)) /
      2.0;
  // RankTwoTensor stress0pos = (_elasticity_tensor[_qp] * projpos * _mechanical_strain[_qp] +
  // (_elasticity_tensor[_qp] * strainpos).initialContraction(IooI))/2.0; RankTwoTensor stress0neg =
  // (_elasticity_tensor[_qp] * projneg * _mechanical_strain[_qp] + (_elasticity_tensor[_qp] *
  // strainneg).initialContraction(IooI))/2.0; RankTwoTensor stress0pos = ((_elasticity_tensor[_qp]
  // * _mechanical_strain[_qp]).initialContraction(projpos) + _elasticity_tensor[_qp] *
  // strainpos)/2.0; RankTwoTensor stress0neg = ((_elasticity_tensor[_qp] *
  // _mechanical_strain[_qp]).initialContraction(projneg) + _elasticity_tensor[_qp] *
  // strainneg)/2.0; RankTwoTensor stress0pos = (_elasticity_tensor[_qp] * projpos *
  // _mechanical_strain[_qp] + _elasticity_tensor[_qp] * strainpos)/2.0; RankTwoTensor stress0neg =
  // (_elasticity_tensor[_qp] * projneg * _mechanical_strain[_qp] + _elasticity_tensor[_qp] *
  // strainneg)/2.0; RankTwoTensor stress0pos =
  // (_mechanical_strain[_qp].initialContraction(_elasticity_tensor[_qp] * projpos) +
  // _elasticity_tensor[_qp] * strainpos)/2.0; RankTwoTensor stress0neg =
  // (_mechanical_strain[_qp].initialContraction(_elasticity_tensor[_qp] * projneg) +
  // _elasticity_tensor[_qp] * strainneg)/2.0; RankTwoTensor stress0pos = ((_elasticity_tensor[_qp]
  // * _mechanical_strain[_qp]).initialContraction(proj_pos) + _elasticity_tensor[_qp] *
  // strainpos)/2.0; RankTwoTensor stress0neg = ((_elasticity_tensor[_qp] *
  // _mechanical_strain[_qp]).initialContraction(proj_neg) + _elasticity_tensor[_qp] *
  // strainneg)/2.0; RankTwoTensor stress0pos = (_elasticity_tensor[_qp] * proj_pos *
  // _mechanical_strain[_qp] + _elasticity_tensor[_qp] * strainpos)/2.0; RankTwoTensor stress0neg =
  // (_elasticity_tensor[_qp] * proj_neg * _mechanical_strain[_qp] + _elasticity_tensor[_qp] *
  // strainneg)/2.0; RankTwoTensor stress0pos = ((_elasticity_tensor[_qp] *
  // proj_pos).inneroutContraction(_mechanical_strain[_qp]) + _elasticity_tensor[_qp] *
  // strainpos)/2.0; RankTwoTensor stress0neg = ((_elasticity_tensor[_qp] *
  // proj_neg).inneroutContraction(_mechanical_strain[_qp]) + _elasticity_tensor[_qp] *
  // strainneg)/2.0;
  // // _stress[_qp] = cfactor * stress0pos * ((1.0 - c) * (1.0 - c) * (1.0 - _kdamage) + _kdamage) + stress0neg;

  // RankTwoTensor stress0neg = (_elasticity_tensor[_qp] * strainneg).initialContraction(proj_neg);
  // RankTwoTensor stress0pos = _elasticity_tensor[_qp] * _mechanical_strain[_qp] - stress0neg;
  // Real G0_neg = (_elasticity_tensor[_qp] * strainneg).doubleContraction(strainneg)/2.0;
  // Real G0_Pos_trial = (_elasticity_tensor[_qp] *
  // _mechanical_strain[_qp]).doubleContraction(_mechanical_strain[_qp])/2.0 - G0_neg;

  Real G0_Pos_trial =
      (_elasticity_tensor[_qp] * strainpos).doubleContraction(_mechanical_strain[_qp]) / 2.0;

  // RankTwoTensor stress0pos = ((_elasticity_tensor[_qp] * IooI).inneroutContraction(strainpos) +
  // (_elasticity_tensor[_qp] * _mechanical_strain[_qp]).chainruleContraction(proj_pos))/2.0;
  // RankTwoTensor stress0neg = ((_elasticity_tensor[_qp] * IooI).inneroutContraction(strainneg) +
  // (_elasticity_tensor[_qp] * _mechanical_strain[_qp]).chainruleContraction(proj_neg))/2.0;

  // RankTwoTensor stress0pos = (strainpos.initialContraction(_elasticity_tensor[_qp] * IooI) +
  // (_elasticity_tensor[_qp] * _mechanical_strain[_qp]).initialContraction(proj_pos))/2.0;
  // RankTwoTensor stress0neg = (strainneg.initialContraction(_elasticity_tensor[_qp] * IooI) +
  // (_elasticity_tensor[_qp] * _mechanical_strain[_qp]).initialContraction(proj_neg))/2.0;
  //
  // Real G0_Pos_trial = (_elasticity_tensor[_qp] *
  // _mechanical_strain[_qp]).doubleContraction(strainpos)/2.0;

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
  _stress[_qp] =
      cfactor * stress0pos * ((1.0 - c) * (1.0 - c) * (1.0 - _kdamage) + _kdamage) + stress0neg;

  // Used in StressDivergencePFFracTensors Jacobian
  _dstress_dc[_qp] = -cfactor * stress0pos * 2.0 * (1.0 - c) * (1.0 - _kdamage);
  // _dstress_dc[_qp] = -cfactor * _dG0_pos_dstrain[_qp] * 2.0 * (1.0 - c);
}

void
AnisoElasticPFDamage::updateJacobian()
{
  Real c = _c[_qp];
  Real cfactor = 1.0;
  RankTwoTensor eigvec;
  RankFourTensor projpos, projneg;
  RankFourTensor naabb, nabab, nabba;
  RankTwoTensor nab, nba, naa, nbb;
  std::vector<Real> eigval(LIBMESH_DIM);
  std::vector<Real> eigabs(LIBMESH_DIM);
  _mechanical_strain[_qp].symmetricEigenvaluesEigenvectors(eigval, eigvec);
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

  RankFourTensor IooI = I.mixedProductIkJl(I);
  // RankFourTensor IoxI = I.outerProduct(I);

  RankTwoTensor eigpos = eigval_tensor * Ipos;
  RankTwoTensor eigneg = eigval_tensor * Ineg;
  // RankTwoTensor strainpos = QoQ * _mechanical_strain[_qp];
  RankTwoTensor strainpos =
      eigvec * eigpos * (eigvec.transpose()); // Postive part of strain == Q*eigpos*QT
  RankTwoTensor strainneg =
      eigvec * eigneg * (eigvec.transpose()); // Negative part of strain == Q*eigpos*QT

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

      Real toldelta = 1e-6; // tolerance value for the eigvalue difference;
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
  RankFourTensor Jpos =
      (_elasticity_tensor[_qp] * projpos * IooI + _elasticity_tensor[_qp] * projpos) / 2.0;
  RankFourTensor Jneg =
      (_elasticity_tensor[_qp] * projneg * IooI + _elasticity_tensor[_qp] * projneg) / 2.0;
  _Jacobian_mult[_qp] = ((1.0 - c) * (1.0 - c) * (1.0 - _kdamage) + _kdamage) * Jpos + Jneg;
  _Jacobian_mult[_qp] = _elasticity_tensor[_qp];
}
