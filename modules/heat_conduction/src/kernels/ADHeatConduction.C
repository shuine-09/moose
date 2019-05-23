//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADHeatConduction.h"

registerMooseObject("HeatConductionApp", ADHeatConduction);

InputParameters
ADHeatConduction::validParams()
{
  InputParameters params = ADDiffusion::validParams();
  params.addParam<MaterialPropertyName>("thermal_conductivity",
                                        "thermal_conductivity",
                                        "the name of the thermal conductivity material property");
  params.addCoupledVar(
      "activated_elem_aux", 1, "Temperature aux variable used to determine activated elements.");
  params.set<bool>("use_displaced_mesh") = true;
  return params;
}

ADHeatConduction::ADHeatConduction(const InputParameters & parameters)
  : ADDiffusion(parameters),
    _thermal_conductivity(getADMaterialProperty<Real>("thermal_conductivity")),
    _activated_elem(coupledValue("activated_elem_aux"))
{
}

ADRealVectorValue
ADHeatConduction::precomputeQpResidual()
{
  if (_activated_elem[_qp] >= 1.0 || _current_elem->subdomain_id() != 1)
    return _thermal_conductivity[_qp] * ADDiffusion::precomputeQpResidual();
  else
    return 0.0;
}
