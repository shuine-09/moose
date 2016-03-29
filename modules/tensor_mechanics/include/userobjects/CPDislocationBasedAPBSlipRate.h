/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CPDISLOCATIONBASEDAPBSLIPRATE_H
#define CPDISLOCATIONBASEDAPBSLIPRATE_H

#include "CrystalPlasticitySlipRate.h"
#include "CrystalPlasticitySlipResistance.h"
#include "RankTwoTensor.h"

class CPDislocationBasedAPBSlipRate;

template<>
InputParameters validParams<CPDislocationBasedAPBSlipRate>();

/**
 * Dislocation based constitutive mode userobject class for apb shear slip rate in the matrix.
 */
class CPDislocationBasedAPBSlipRate : public CrystalPlasticitySlipRate
{
 public:
  CPDislocationBasedAPBSlipRate(const InputParameters & parameters);

  virtual bool calcSlipRate(unsigned int qp, Real dt, std::vector<Real> & val) const;
  virtual bool calcSlipRateDerivative(unsigned int qp, Real dt, std::vector<Real> & val) const;
  virtual void calcFlowDirection(unsigned int qp, std::vector<RankTwoTensor> & flow_direction) const;

 protected:
  const MaterialProperty<std::vector<Real> > & _mat_prop_mobile_dislocation_density;

  const MaterialProperty<std::vector<Real> > & _mat_prop_immobile_dislocation_density;

  const MaterialProperty<std::vector<Real> > & _mat_prop_thermal_slip_resistance;

  const MaterialProperty<std::vector<Real> > & _mat_prop_apb_slip_resistance;

  const MaterialProperty<RankTwoTensor> & _pk2;

  const MaterialProperty<std::vector<RankTwoTensor> > & _flow_direction;

  Real _b, _lg, _jump_freq, _enthal;
  Real _k, _temp;
  Real _factor;
  Real _critial_density;

};

#endif // CPDISLOCATIONBASEDAPBSLIPRATE_H
