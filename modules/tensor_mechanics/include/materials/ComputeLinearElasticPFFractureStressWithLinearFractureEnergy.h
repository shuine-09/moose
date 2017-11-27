//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef COMPUTELINEARELASTICPFFRACTURESTRESSWITHLINEARFRACTUREENERGY_H
#define COMPUTELINEARELASTICPFFRACTURESTRESSWITHLINEARFRACTUREENERGY_H

#include "ComputeIsotropicLinearElasticPFFractureStress.h"

class ComputeLinearElasticPFFractureStressWithLinearFractureEnergy;

template <>
InputParameters validParams<ComputeLinearElasticPFFractureStressWithLinearFractureEnergy>();

/**
 * Phase-field fracture
 * This class computes the stress and energy contribution to fracture
 * Small strain Anisotropic Elastic formulation
 * Stiffness matrix scaled for heterogeneous elasticity property
 */
class ComputeLinearElasticPFFractureStressWithLinearFractureEnergy
  : public ComputeIsotropicLinearElasticPFFractureStress
{
public:
  ComputeLinearElasticPFFractureStressWithLinearFractureEnergy(const InputParameters & parameters);

protected:
  virtual void computeQpStress();
  /// Function required to initialize statefull material properties
  virtual void initQpStatefulProperties();
};

#endif // COMPUTELINEARELASTICPFFRACTURESTRESSWITHLINEARFRACTUREENERGY_H
