//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef COMPUTEWELDTHERMALEXPANSIONEIGENSTRAIN_H
#define COMPUTEWELDTHERMALEXPANSIONEIGENSTRAIN_H

#include "Material.h"
#include "ComputeThermalExpansionEigenstrainBase.h"
#include "DerivativeMaterialInterface.h"

class ComputeWeldThermalExpansionEigenstrain;
class WeldStateIndicator;
class Function;

template <>
InputParameters validParams<ComputeWeldThermalExpansionEigenstrain>();

/**
 * ComputeWeldThermalExpansionEigenstrain computes an eigenstrain for thermal expansion
 * for after the material has solidified using a constant expansion coefficient.
 */
class ComputeWeldThermalExpansionEigenstrain : public ComputeThermalExpansionEigenstrainBase
{
public:
  ComputeWeldThermalExpansionEigenstrain(const InputParameters & parameters);

protected:
  virtual void computeThermalStrain(Real & thermal_strain, Real & instantaneous_cte) override;

  const Real & _thermal_expansion_coeff;
  const WeldStateIndicator * _weld_state;
};

#endif // COMPUTEWELDTHERMALEXPANSIONEIGENSTRAIN_H
