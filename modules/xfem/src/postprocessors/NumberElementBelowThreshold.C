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

#include "NumberElementBelowThreshold.h"

template <>
InputParameters
validParams<NumberElementBelowThreshold>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  return params;
}

NumberElementBelowThreshold::NumberElementBelowThreshold(const InputParameters & parameters)
  : ElementIntegralPostprocessor(parameters), _weibull(getMaterialProperty<Real>("weibull"))
{
}

Real
NumberElementBelowThreshold::computeIntegral()
{
  if (_weibull[0] <= 100)
    return 1;
  else
    return 0;
}

Real
NumberElementBelowThreshold::computeQpIntegral()
{
  return 0;
}
