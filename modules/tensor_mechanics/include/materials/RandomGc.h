//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#ifndef RANDOMMGC_H
#define RANDOMMGC_H

#include "Material.h"

class RandomGc;

template <>
InputParameters validParams<RandomGc>();

class RandomGc : public Material
{
public:
  RandomGc(const InputParameters & parameters);
  virtual void computeQpProperties() override;
  virtual void initQpStatefulProperties() override;

protected:
  /// Material property defining gc parameter, declared elsewhere
  MaterialProperty<Real> & _gc;
  const MaterialProperty<Real> & _gc_old;

  Real _gc_mean;
  Real _random_range;
};

#endif // RANDOMMGC_H
