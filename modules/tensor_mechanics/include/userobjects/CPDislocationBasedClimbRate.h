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
#include "Function.h"

class CPDislocationBasedClimbRate;

template <>
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
  virtual void calcFlowDirection(unsigned int qp,
                                 std::vector<RankTwoTensor> & flow_direction) const;

protected:
  // const MaterialProperty<std::vector<Real>> & _cv;
  const MaterialProperty<std::vector<Real>> & _rho_m;
  const MaterialProperty<std::vector<Real>> & _rho_i;
  const MaterialProperty<RankTwoTensor> & _pk2;
  const MaterialProperty<std::vector<RankTwoTensor>> & _flow_direction;

  Real _diffusivity;
  Real _rc;
  Real _b;
  Real _molar_volume;
  Real _atom_volume;
  Real _lattice_constant;
  Real _boltz_const;
  Real _gas_constant;
  Real _temp;
  Real _stress_factor;
  Real _diffusivity_factor;
  Real _activation_energy;
  Real _precipitate_radius;
  Real _precipitate_volume_fraction;
  Real _theta;
  Real _use_damage;
  const MaterialProperty<Real> * _w_old;
};

#endif // CPDISLOCATIONBASEDCLIMB_H
