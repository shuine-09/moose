/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CRYSTALPLASTICITYSLIPRATEGSS_H
#define CRYSTALPLASTICITYSLIPRATEGSS_H

#include "CrystalPlasticitySlipRate.h"
#include "RankTwoTensor.h"

class CrystalPlasticitySlipRateGSS;

template<>
InputParameters validParams<CrystalPlasticitySlipRateGSS>();

/**
 * phenomenological constitutive models' slip rate userobject class
 * The virtual functions written below must be
 * over-ridden in derived classes to provide actual values
 */

class CrystalPlasticitySlipRateGSS : public CrystalPlasticitySlipRate
{
 public:
  CrystalPlasticitySlipRateGSS(const InputParameters & parameters);

  virtual bool calcSlipRate(unsigned qp, Real dt, std::vector<RankTwoTensor> & schmid_tensor, std::vector<Real> & val) const;

  virtual bool calcSlipRateDerivative(unsigned qp, Real dt, std::vector<RankTwoTensor> & schmid_tensor, std::vector<Real> & val) const;

  virtual void calcSchmidTensor(RankTwoTensor & crysrot, std::vector<RankTwoTensor> & schmid_tensor) const;

 protected:
  
  virtual void readFileFlowRateParams();

  virtual void getFlowRateParams();

  std::vector<const MaterialProperty<std::vector<Real> > * > _mat_prop_state_var_local;
  
  unsigned int _num_mat_state_var;

  const MaterialProperty<RankTwoTensor> & _pk2;
  
  DenseVector<Real> _a0;
  DenseVector<Real> _xm;
};

#endif // CRYSTALPLASTICITYSLIPRATEGSS_H
