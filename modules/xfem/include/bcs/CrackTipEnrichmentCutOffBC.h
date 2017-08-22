/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef CRACKTIPENRICHMENTCUTOFFBC_H
#define CRACKTIPENRICHMENTCUTOFFBC_H

#include "PresetBC.h"
#include "CrackFrontDefinition.h"

class CrackTipEnrichmentCutOffBC;

template <>
InputParameters validParams<CrackTipEnrichmentCutOffBC>();

class CrackTipEnrichmentCutOffBC : public PresetBC
{
public:
  CrackTipEnrichmentCutOffBC(const InputParameters & parameters);

protected:
  virtual bool shouldApply() override;

  const Real _cut_off_radius;

private:
  const CrackFrontDefinition * _crack_front_definition;
};

#endif /* CRACKTIPENRICHMENTCUTOFFBC_H */
