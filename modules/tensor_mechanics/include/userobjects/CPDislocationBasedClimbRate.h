/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CPDISLOCATIONBASEDCLIMBRATE_H
#define CPDISLOCATIONBASEDCLIMBRATE_H

#include "CrystalPlasticitySlipRate.h"
#include "CrystalPlasticitySlipResistance.h"
#include "RankTwoTensor.h"

class CPDislocationBasedClimbRate;

template<>
InputParameters validParams<CPDislocationBasedClimbRate>();

/**
 * Dislocation based constitutive mode userobject class for climb rate.
 */
class CPDislocationBasedClimbRate : public CrystalPlasticitySlipRate
{
 public:
  CPDislocationBasedClimbRate(const InputParameters & parameters);

  virtual bool calcSlipRate(unsigned int qp, Real dt, std::vector<Real> & val) const;
  virtual bool calcSlipRateDerivative(unsigned int qp, Real dt, std::vector<Real> & val) const;
  virtual void calcFlowDirection(unsigned int qp, std::vector<RankTwoTensor> & flow_direction) const;

 protected:

  const MaterialProperty<std::vector<Real> > & _rho_m;
  const MaterialProperty<std::vector<Real> > & _rho_i;
  const MaterialProperty<RankTwoTensor> & _pk2;
  const MaterialProperty<std::vector<RankTwoTensor> > & _flow_direction;

  Real _climb_resistance;
  Real _b;
  Real _l;
  Real _p;
  Real _q;
  Real _enthal;
  Real _k;
  Real _temp;
  Real _edge_screw_angle;
};

#endif // CPDISLOCATIONBASEDCLIMB_H
