//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef COMPUTEWELDINGISOTROPICELASTICITYTENSOR_H
#define COMPUTEWELDINGISOTROPICELASTICITYTENSOR_H

#include "ComputeIsotropicElasticityTensor.h"

class ComputeWeldingIsotropicElasticityTensor;
class WeldStateIndicator;

template <>
InputParameters validParams<ComputeIsotropicElasticityTensor>();

/**
 * ComputeIsotropicElasticityTensor defines an elasticity tensor material for
 * isotropic materials.
 */
class ComputeWeldingIsotropicElasticityTensor : public ComputeIsotropicElasticityTensor
{
public:
  ComputeWeldingIsotropicElasticityTensor(const InputParameters & parameters);

protected:
  virtual void computeQpElasticityTensor() override;

  const WeldStateIndicator * _weld_state;
};

#endif // COMPUTEWELDINGISOTROPICELASTICITYTENSOR_H
