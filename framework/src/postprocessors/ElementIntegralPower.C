/*************************************************/
/*           DO NOT MODIFY THIS HEADER           */
/*                                               */
/*                     BISON                     */
/*                                               */
/*    (c) 2015 Battelle Energy Alliance, LLC     */
/*            ALL RIGHTS RESERVED                */
/*                                               */
/*   Prepared by Battelle Energy Alliance, LLC   */
/*     Under Contract No. DE-AC07-05ID14517      */
/*     With the U. S. Department of Energy       */
/*                                               */
/*     See COPYRIGHT for full restrictions       */
/*************************************************/

#include "ElementIntegralPower.h"

#include "BurnupFunction.h"

template<>
InputParameters validParams<ElementIntegralPower>()
{
  InputParameters params = validParams<ElementIntegralVariablePostprocessor>();
  params.addCoupledVar("fission_rate", "Coupled fission rate");
  params.addParam<FunctionName>("burnup_function", "Burnup function");
  params.addParam<Real>("energy_per_fission", 3.28451e-11, "Energy released per fission (J/fission)");
  return params;
}

ElementIntegralPower::ElementIntegralPower(const InputParameters & parameters) :
    ElementIntegralVariablePostprocessor(parameters),
    _energy_per_fission(getParam<Real>("energy_per_fission")),
    _has_fission_rate(isCoupled("fission_rate")),
    _fission_rate(_has_fission_rate?coupledValue("fission_rate"):_zero),
    _burnup_function( isParamValid("burnup_function") ?
                      dynamic_cast<BurnupFunction*>(&getFunction("burnup_function")) : NULL )

{
  if (_has_fission_rate && _burnup_function)
  {
    mooseError("Cannot have fission_rate and burnup_function in ElementIntegralPower");
  }
}

Real
ElementIntegralPower::computeQpIntegral()
{
  const Real fission_rate = _has_fission_rate ? _fission_rate[_qp] : _burnup_function->fissionRate(_q_point[_qp]);
  return fission_rate * _energy_per_fission;
}
