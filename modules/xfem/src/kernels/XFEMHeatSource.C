/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "XFEMHeatSource.h"
// MOOSE
#include "Function.h"


template<>
InputParameters validParams<XFEMHeatSource>()
{
  InputParameters params = validParams<BodyForce>();

  // Override defaults and documentation, weak form is identical to BodyForce in MOOSE
  params.addParam<Real>("value", 1.0, "Value of heat source. Multiplied by function if present.");
  params.addParam<FunctionName>("function", "1", "Function describing the volumetric heat source");
  return params;
}

XFEMHeatSource::XFEMHeatSource(const InputParameters & parameters) :
    BodyForce(parameters)
{
}

Real
XFEMHeatSource::computeQpResidual()
{
  Real factor = _value * _function.value(_t, _q_point[_qp]);

  if (_q_point[_qp](1) > 1.0)
    factor = 0.0;

  if (_postprocessor)
    factor *= *_postprocessor;
  return _test[_i][_qp] * -factor;
}
