/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CPDISLOCATIONBASEDGLIDESLIPRATE_H
#define CPDISLOCATIONBASEDGLIDESLIPRATE_H

#include "CrystalPlasticitySlipRate.h"
#include "CrystalPlasticitySlipResistance.h"
#include "RankTwoTensor.h"

class CPDislocationBasedGlideSlipRate;

template<>
InputParameters validParams<CPDislocationBasedGlideSlipRate>();

/**
 * Dislocation based constitutive mode userobject class for glide slip rate in the matrix.
 */
class CPDislocationBasedGlideSlipRate : public CrystalPlasticitySlipRate
{
 public:
  CPDislocationBasedGlideSlipRate(const InputParameters & parameters);

  virtual bool calcSlipRate(unsigned int qp, Real dt, std::vector<Real> & val) const;
  virtual bool calcSlipRateDerivative(unsigned int qp, Real dt, std::vector<Real> & val) const;
  virtual void calcFlowDirection(unsigned int qp, std::vector<RankTwoTensor> & flow_direction) const;

 protected:
  const MaterialProperty<std::vector<Real> > & _mat_prop_mobile_dislocation_density;

  const MaterialProperty<std::vector<Real> > & _mat_prop_thermal_slip_resistance;

  const MaterialProperty<std::vector<Real> > & _mat_prop_athermal_slip_resistance;

  const MaterialProperty<RankTwoTensor> & _pk2;

  const MaterialProperty<std::vector<RankTwoTensor> > & _flow_direction;

  ///Penalty parameter value used to regularize activation energy based flow rule.
  Real _penalty_param;

  Real _b, _lg, _jump_freq, _enthal, _p, _q;
  Real _k, _temp;

};

#endif // CPDISLOCATIONBASEDGLIDESLIPRATE_H
