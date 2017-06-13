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
    _G0_pos_old(declarePropertyOld<Real>("G0_pos")),  //added G0_old variable for old time step
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
  if (_t > 0)
  _G0_pos_old[_qp] = _G0_pos[_qp];
  else{
    _G0_pos[_qp] = 0.0;
    _G0_pos_old[_qp] = 0.0;
  }
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
  if (c > 1.0)
    cfactor = 0.0;

  _Jacobian_mult[_qp] = _elasticity_tensor[_qp];

  RankTwoTensor eigvec;
  std::vector<Real> eigval(LIBMESH_DIM);
  _mechanical_strain[_qp].symmetricEigenvaluesEigenvectors(eigval, eigvec);
  RankTwoTensor Ipos, Ineg, eigval_tensor;
  RankTwoTensor I(RankTwoTensor::initIdentity);
  eigval_tensor.fillFromInputVector(eigval);
  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    if (eigval[i] < 0.0)
      Ineg(i, i) = 1.0;

  Ipos = I - Ineg;
  // Creation of some critical fourth order tensors
  RankFourTensor QoQ = eigvec.mixedProductIkJl(eigvec);
  RankFourTensor QToQT = (eigvec.transpose()).mixedProductIkJl(eigvec.transpose());
  RankFourTensor IpoIp = Ipos.mixedProductIkJl(I);
  RankFourTensor InoIn = Ineg.mixedProductIkJl(I);
  // RankFourTensor QIoQI = (eigvec.transpose()).outerProduct(eigvec.transpose());


  // const RankFourTensor Jnc = _Jacobian_mult[_qp] * QoQ * IpoIp * QToQT; //Jacobian matrix under new coordinate system
  const RankFourTensor Jnc = QToQT * _elasticity_tensor[_qp] * QoQ; //Jacobian matrix under new coordinate system Not sure
  // const RankFourTensor Jnc = QoQ * _Jacobian_mult[_qp] * QToQT; //Jacobian matrix under new coordinate system Not sure
  // const RankFourTensor Jnc = _Jacobian_mult[_qp] * QoQ * QToQT; //Jacobian matrix under new coordinate system  //Not sure this is correct


  const RankFourTensor Jpos = _elasticity_tensor[_qp] * QoQ * IpoIp * QToQT;
  const RankFourTensor Jneg = _elasticity_tensor[_qp] * QoQ * InoIn * QToQT;
  _Jacobian_mult[_qp] = cfactor * Jpos * ((1.0 - c) * (1.0 - c) + _kdamage) + Jneg;

  RankTwoTensor eigpos = eigval_tensor * Ipos;  //positive principle strain tensor
  RankTwoTensor eigneg = eigval_tensor * Ineg;  //negative eigen strain tensor

  RankTwoTensor stressncpos = Jnc * eigpos;
  RankTwoTensor stressncneg = Jnc * eigneg;

  // RankTwoTensor stressncmark = Jnc * eigval_tensor;
  // RankTwoTensor stressncneg = Jnc * eigneg;

  RankTwoTensor stress0pos = eigvec * stressncpos * (eigvec.transpose());
  RankTwoTensor stress0neg = eigvec * stressncneg * (eigvec.transpose());

  // RankTwoTensor stress0pos = QIoQI * stressncpos;
  // RankTwoTensor stress0neg = QIoQI * stressncneg;

  // Damage associated with positive component of stress
  _stress[_qp] = cfactor * stress0pos * ((1.0 - c) * (1.0 - c) + _kdamage) + stress0neg;

  // Energy with positive principal strains
  _G0_pos[_qp] = stressncpos.doubleContraction(eigval_tensor) / 2.0;
  // _G0_neg[_qp] = stressncneg.doubleContraction(eigval_tensor) / 2.0;



  //Before Update
  // Real G0_Pos_trial = stressncmark.doubleContraction(eigval_tensor) / 2.0;
  Real G0_Pos_trial = stressncpos.doubleContraction(eigval_tensor) / 2.0;
  // Real G0_Pos_trial = stressncpos.doubleContraction(eigpos) / 2.0;


  if (G0_Pos_trial > _G0_pos_old[_qp]){
    _G0_pos[_qp] = G0_Pos_trial;
    _dG0_pos_dstrain[_qp] = stress0pos;
  }else{
    _G0_pos[_qp] = _G0_pos_old[_qp];
    _dG0_pos_dstrain[_qp] = stress0pos * 0.0;  // derivative of positive strain energy wrt strain
  }

  // Used in PFFracBulkRate Jacobian
  // _dG0_pos_dstrain[_qp] = stress0pos;

  // Damage associated with positive component of stress
  _stress[_qp] = cfactor * stress0pos * ((1.0 - c) * (1.0 - c) + _kdamage) + stress0neg;

  // Used in StressDivergencePFFracTensors Jacobian
  _dstress_dc[_qp] = -cfactor * stress0pos * 2.0 * (1.0 - c);
  // _dstress_dc[_qp] = -cfactor * _dG0_pos_dstrain[_qp] * 2.0 * (1.0 - c);
}

void
AnisoElasticPFDamage::updateJacobian()
{
  _Jacobian_mult[_qp] = _elasticity_tensor[_qp];
}
