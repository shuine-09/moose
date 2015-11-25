/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CRYSTALPLASTICITYSLIPRESISTANCEGSS_H
#define CRYSTALPLASTICITYSLIPRESISTANCEGSS_H

#include "CrystalPlasticitySlipResistance.h"
#include "RankTwoTensor.h"

class CrystalPlasticitySlipResistanceGSS;

template<>
InputParameters validParams<CrystalPlasticitySlipResistanceGSS>();

/**
 * phenomenological constitutive models' slip resistance userobject class
 * The virtual functions written below must be
 * over-ridden in derived classes to provide actual values
 */

class CrystalPlasticitySlipResistanceGSS : public CrystalPlasticitySlipResistance
{
 public:

  CrystalPlasticitySlipResistanceGSS(const InputParameters & parameters);

  virtual bool calcSlipResistance(unsigned int qp, std::vector<Real> & val) const;
 
  virtual void initSlipSysProps(std::vector<Real> & val) const;

  virtual void assignSlipSysRes(std::vector<Real> & val) const;

  virtual void readFileInitSlipSysRes(std::vector<Real> & val) const;

  virtual void getInitSlipSysRes(std::vector<Real> & val) const;

  virtual std::vector<Real> getHardnessParams() const;

 protected:

  virtual void readFileHardnessParams();

  virtual void assignHardnessParams();

  ///File should contain initial values of the slip system resistances.
  std::string _slip_sys_res_prop_file_name;

  ///The hardening parameters in this class are read from .i file. The user can override to read from file.
  std::string _slip_sys_hard_prop_file_name;

  ///Read from options for initial values of internal variables
  MooseEnum _intvar_read_type;

  std::vector<Real> _gprops;
  std::vector<Real> _hprops;

  std::vector<const MaterialProperty<std::vector<Real> > * > _mat_prop_state_var_local;
  unsigned int _num_mat_state_var;

  const MaterialProperty<std::vector<Real> > & _mat_prop_slip_rate;

  Real _h0;
  Real _tau_sat;
  Real _tau_init;
  Real _r;
};

#endif //  CRYSTALPLASTICITYSLIPRESISTANCEGSS_H
