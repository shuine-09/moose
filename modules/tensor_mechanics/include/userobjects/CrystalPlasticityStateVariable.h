/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CRYSTALPLASTICITYSTATEVARIABLE_H
#define CRYSTALPLASTICITYSTATEVARIABLE_H

#include "CrystalPlasticityUOBase.h"
#include "RankTwoTensor.h"

class CrystalPlasticityStateVariable;

template<>
InputParameters validParams<CrystalPlasticityStateVariable>();

/**
 * Crystal plasticity state variable userobject class
 * The virtual functions written below must be
 * over-ridden in derived classes to provide actual values
 */

class CrystalPlasticityStateVariable : public CrystalPlasticityUOBase
{
 public:
  CrystalPlasticityStateVariable(const InputParameters & parameters);

  /// Returns the slip rate  name  (eg "DislocationGlide")
  virtual std::string crystalPlasticityUOName() const;

  virtual bool updateStateVariable(unsigned int qp, Real dt, std::vector<Real> & val) const = 0;

  virtual void initSlipSysProps(std::vector<Real> & val) const = 0;

  virtual void assignSlipSysRes(std::vector<Real> & val){};

  virtual void readFileInitSlipSysRes(std::vector<Real> & val){};

  virtual void getInitSlipSysRes(std::vector<Real> & val){};

 protected:

  unsigned int _nss;
};

#endif // CRYSTALPLASTICITYSTATEVARIABLE_H
