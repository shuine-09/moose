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

#include "PresetBCNearTipEnrichment.h"

template<>
InputParameters validParams<PresetBCNearTipEnrichment>()
{
  InputParameters p = validParams<NodalBC>();
  p.addRequiredParam<Real>("value", "Value of the BC");
  p.addParam<Real>("radius", 0.1, "Radius");
  return p;
}


PresetBCNearTipEnrichment::PresetBCNearTipEnrichment(const InputParameters & parameters) :
  PresetNodalBC(parameters),
  _value(getParam<Real>("value")),
  _radius(getParam<Real>("radius"))
{
}

Real
PresetBCNearTipEnrichment::computeQpValue()
{
  return _value;
}

bool
PresetBCNearTipEnrichment::shouldApply()
{
  Real dist = std::sqrt(((*_current_node)(0) - 0.5) * ((*_current_node)(0) - 0.5) + ((*_current_node)(1) - 0.5) * ((*_current_node)(1) - 0.5));

  if (dist > _radius)
    return true;
  else
    return false;
//  return true;
}
