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

#ifndef PRESETBCNEARTIPENRICHMENT_H
#define PRESETBCNEARTIPENRICHMENT_H

#include "PresetNodalBC.h"

class PresetBCNearTipEnrichment;

template<>
InputParameters validParams<PresetBCNearTipEnrichment>();

class PresetBCNearTipEnrichment : public PresetNodalBC
{
public:
  PresetBCNearTipEnrichment(const InputParameters & parameters);

protected:
  virtual Real computeQpValue() override;

  virtual bool shouldApply() override;

  const Real & _value;
  const Real _radius;
};

#endif /* PRESETBCNEARTIPENRICHMENT_H */
