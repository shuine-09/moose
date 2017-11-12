/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef VESICLEINITIALAREAVOLUMEUO_H
#define VESICLEINITIALAREAVOLUMEUO_H

#include "DiscreteElementUserObject.h"

class VesicleInitialAreaVolumeUO;

template<>
InputParameters validParams<VesicleInitialAreaVolumeUO>();

class VesicleInitialAreaVolumeUO : public DiscreteElementUserObject
{
 public:
  VesicleInitialAreaVolumeUO(const InputParameters & parameters);

  virtual void setInitialAreaVolume(Real area, Real volume);

  virtual Real initialArea() const;

  virtual Real initialVolume() const;

 protected:
  Real _volume0;
  Real _area0;
};

#endif // VESICLEINITIALAREAVOLUMEUO_H
