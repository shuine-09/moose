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

#include "VesicleAreaPostprocessor.h"

registerMooseObject("MooseApp", VesicleAreaPostprocessor);

template <>
InputParameters
validParams<VesicleAreaPostprocessor>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredCoupledVar("variable",
                               "The name of the variable that this userobject applies to");
  params.addParam<Real>("epsilon", 0.01, "The interfacial penalty parameter.");
  return params;
}

VesicleAreaPostprocessor::VesicleAreaPostprocessor(const InputParameters & parameters)
  : ElementIntegralPostprocessor(parameters),
    MooseVariableInterface<Real>(this, false),
    _u(coupledValue("variable")),
    _grad_u(coupledGradient("variable")),
    _epsilon(getParam<Real>("epsilon")),
    _qp(0)
{
}

void
VesicleAreaPostprocessor::threadJoin(const UserObject & y)
{
  const VesicleAreaPostprocessor & pps = static_cast<const VesicleAreaPostprocessor &>(y);
  _integral_value += pps._integral_value;
}

Real
VesicleAreaPostprocessor::computeQpIntegral()
{
  Real value = 0.5 * _epsilon * _grad_u[_qp] * _grad_u[_qp] +
               1.0 / (4.0 * _epsilon) * pow((_u[_qp] * _u[_qp] - 1.0), 2.0);
  // Real value = _epsilon * _grad_u[_qp] * _grad_u[_qp];
  value *= 3.0 / 2.0 / std::sqrt(2.0);
  return value;
}
