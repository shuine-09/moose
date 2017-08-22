/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "CrackTipEnrichmentCutOffBC.h"

template <>
InputParameters
validParams<CrackTipEnrichmentCutOffBC>()
{
  InputParameters p = validParams<PresetBC>();
  p.addParam<Real>("cut_off_radius", 0.1, "Radius");
  p.set<bool>("use_displaced_mesh") = false;
  return p;
}

CrackTipEnrichmentCutOffBC::CrackTipEnrichmentCutOffBC(const InputParameters & parameters)
  : PresetBC(parameters), _cut_off_radius(getParam<Real>("cut_off_radius"))
{
}

bool
CrackTipEnrichmentCutOffBC::shouldApply()
{
  // TODO: crack tip coordinate
  Real dist = std::sqrt(((*_current_node)(0) - 0.5) * ((*_current_node)(0) - 0.5) +
                        ((*_current_node)(1) - 1.0) * ((*_current_node)(1) - 1.0));

  if (dist > _cut_off_radius)
    return true;
  else
    return false;
}
