/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "VesicleInitialAreaVolumeUO.h"

registerMooseObject("MooseApp", VesicleInitialAreaVolumeUO);

template <>
InputParameters
validParams<VesicleInitialAreaVolumeUO>()
{
  InputParameters params = validParams<DiscreteElementUserObject>();
  return params;
}

VesicleInitialAreaVolumeUO::VesicleInitialAreaVolumeUO(const InputParameters & parameters)
  : DiscreteElementUserObject(parameters)
{
}

void
VesicleInitialAreaVolumeUO::setInitialAreaVolume(Real area, Real volume)
{
  _area0 = area;
  _volume0 = volume;
}

Real
VesicleInitialAreaVolumeUO::initialArea() const
{
  return _area0;
}

Real
VesicleInitialAreaVolumeUO::initialVolume() const
{
  return _volume0;
}
