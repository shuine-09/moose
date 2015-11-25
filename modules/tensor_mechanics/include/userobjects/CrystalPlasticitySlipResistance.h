/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CRYSTALPLASTICITYSLIPRESISTANCE_H
#define CRYSTALPLASTICITYSLIPRESISTANCE_H

#include "CrystalPlasticityUOBase.h"
#include "RankTwoTensor.h"

class CrystalPlasticitySlipResistance;

template<>
InputParameters validParams<CrystalPlasticitySlipResistance>();

/**
 * Crystal plasticity slip resistance userobject class
 * The virtual functions written below must be
 * over-ridden in derived classes to provide actual values
 */

class CrystalPlasticitySlipResistance : public CrystalPlasticityUOBase
{
 public:
  CrystalPlasticitySlipResistance(const InputParameters & parameters);

  virtual std::string crystalPlasticityUOName() const;
  
  virtual bool calcSlipResistance(unsigned int qp, std::vector<Real> & val) const = 0;

  virtual void initSlipSysProps(std::vector<Real> & val) const;

  virtual std::vector<Real> getHardnessParams() const;

 protected:

  unsigned int _nss;
};

#endif // CRYSTALPLASTICITYSLIPRESISTANCE_H
