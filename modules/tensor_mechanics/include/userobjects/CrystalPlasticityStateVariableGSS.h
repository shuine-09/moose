/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CRYSTALPLASTICITYSTATEVARIABLEGSS_H
#define CRYSTALPLASTICITYSTATEVARIABLEGSS_H

#include "CrystalPlasticityStateVariable.h"
#include "RankTwoTensor.h"
#include "CrystalPlasticitySlipResistance.h"

class CrystalPlasticityStateVariableGSS;

template<>
InputParameters validParams<CrystalPlasticityStateVariableGSS>();

/**
 * Crystal plasticity state variable userobject class
 * The virtual functions written below must be
 * over-ridden in derived classes to provide actual values
 */

class CrystalPlasticityStateVariableGSS : public CrystalPlasticityStateVariable 
{
 public:

   CrystalPlasticityStateVariableGSS(const InputParameters & parameters);

   virtual bool updateStateVariable(unsigned int qp, Real dt, std::vector<Real> & val) const;

   virtual void initSlipSysProps(std::vector<Real> & val) const;

 protected:

   bool _use_slip_resistance_as_state_var;

   const CrystalPlasticitySlipResistance * _uo_slip_resistance;

   unsigned int _num_mat_state_var_evol_rates;

   unsigned int _num_mat_state_var;
   
   std::vector<const MaterialProperty<std::vector<Real> > * > _mat_prop_state_var_evol_rates;

   std::vector<const MaterialProperty<std::vector<Real> > * > _mat_prop_state_var_local;

   std::vector<const MaterialProperty<std::vector<Real> > * > _mat_prop_state_var_local_old;
};

#endif // CRYSTALPLASTICITYSTATEVARIABLE_H
