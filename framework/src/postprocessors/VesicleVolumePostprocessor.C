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

#include "VesicleVolumePostprocessor.h"

registerMooseObject("MooseApp", VesicleVolumePostprocessor);

template <>
InputParameters
validParams<VesicleVolumePostprocessor>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredCoupledVar("variable",
                               "The name of the variable that this userobject applies to");
  return params;
}

VesicleVolumePostprocessor::VesicleVolumePostprocessor(const InputParameters & parameters)
  : ElementIntegralPostprocessor(parameters),
    MooseVariableInterface<Real>(this, false),
    _u(coupledValue("variable")),
    _grad_u(coupledGradient("variable")),
    _qp(0)
{
}

void
VesicleVolumePostprocessor::threadJoin(const UserObject & y)
{
  const VesicleVolumePostprocessor & pps = static_cast<const VesicleVolumePostprocessor &>(y);
  _integral_value += pps._integral_value;
}

Real
VesicleVolumePostprocessor::computeQpIntegral()
{
  return 0.5 * (1.0 - _u[_qp]);
}
