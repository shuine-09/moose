/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CPDISLOCATIONBASEDCLIMBRATEGEERS_H
#define CPDISLOCATIONBASEDCLIMBRATEGEERS_H

#include "CrystalPlasticitySlipRate.h"
#include "CrystalPlasticitySlipResistance.h"
#include "RankTwoTensor.h"

class CPDislocationBasedClimbRateGeers;

template <>
InputParameters validParams<CPDislocationBasedClimbRateGeers>();

/**
 * Dislocation based constitutive mode userobject class for climb rate.
 */
class CPDislocationBasedClimbRateGeers : public CrystalPlasticitySlipRate
{
public:
  CPDislocationBasedClimbRateGeers(const InputParameters & parameters);

  virtual bool calcSlipRate(unsigned int qp, Real dt, std::vector<Real> & val) const;
  virtual bool calcSlipRateDerivative(unsigned int qp, Real dt, std::vector<Real> & val) const;
  virtual void calcFlowDirection(unsigned int qp,
                                 std::vector<RankTwoTensor> & flow_direction) const;

protected:
  const MaterialProperty<std::vector<Real>> & _cv;
  const MaterialProperty<std::vector<Real>> & _rho_m;
  const MaterialProperty<std::vector<Real>> & _rho_i;
  const MaterialProperty<RankTwoTensor> & _pk2;
  const MaterialProperty<std::vector<RankTwoTensor>> & _flow_direction;

  Real _c0;
  Real _diffusivity;
  Real _rc;
  Real _b;
  Real _molar_volume;
  Real _gas_constant;
  Real _temp;
};

#endif // CPDISLOCATIONBASEDCLIMBGEERS_H
