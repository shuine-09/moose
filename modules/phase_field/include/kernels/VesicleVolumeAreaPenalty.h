/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef VESICLEVOLUMEAREAPENALTY_H
#define VESICLEVOLUMEAREAPENALTY_H

#include "Kernel.h"
#include "VesicleVolumePostprocessor.h"
#include "VesicleAreaPostprocessor.h"
#include "VesicleInitialAreaVolumeUO.h"

class VesicleVolumeAreaPenalty;

template <>
InputParameters validParams<VesicleVolumeAreaPenalty>();

class VesicleVolumeAreaPenalty : public Kernel
{
public:
  VesicleVolumeAreaPenalty(const InputParameters & parameters);

  virtual ~VesicleVolumeAreaPenalty();

protected:
  virtual void timestepSetup();

  virtual void jacobianSetup();

  virtual void residualSetup();

  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  Real _alpha_v;

  Real _alpha_a;

  Real _epsilon;
  bool _rz;

  Real _volume, _volume_0;
  Real _area, _area_0;

  const VariablePhiSecond & _second_phi;
  const VariableTestSecond & _second_test;
  const VariableSecond & _second_u;

  const PostprocessorValue & _vesicle_area;
  const PostprocessorValue & _vesicle_volume;

  bool _use_prescribed_volume;
  bool _use_prescribed_area;

  Real _prescribed_volume;
  Real _prescribed_area;

  Real _alpha_v0;
  Real _alpha_a0;

  bool _use_nonlocal_constraint;

  // const VesicleInitialAreaVolumeUO & _uo_initial_area_volume;
};

#endif /* VESICLEVOLUMEAREAPENALTY_H */
