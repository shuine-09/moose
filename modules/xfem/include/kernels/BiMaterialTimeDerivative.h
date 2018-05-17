//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef BIMATERIALTIMEDERIVATIVE_H
#define BIMATERIALTIMEDERIVATIVE_H

#include "TimeDerivative.h"
#include "Material.h"

// Forward Declaration
class BiMaterialTimeDerivative;

template <>
InputParameters validParams<BiMaterialTimeDerivative>();

class BiMaterialTimeDerivative : public TimeDerivative
{
public:
  BiMaterialTimeDerivative(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  const MaterialProperty<Real> & _time_step_scale;
};

#endif // BIMATERIALTIMEDERIVATIVE_H
