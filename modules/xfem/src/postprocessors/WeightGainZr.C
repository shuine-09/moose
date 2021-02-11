//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "WeightGainZr.h"

registerMooseObject("MooseApp", WeightGainZr);

defineLegacyParams(WeightGainZr);

InputParameters
WeightGainZr::validParams()
{
  InputParameters params = GeneralPostprocessor::validParams();
  params.addParam<Real>("temperature",1473.15,"Temperature of the cladding (K)");
  params.addParam<PostprocessorName>("flux_integral", "The name of the vacancy flux inegral postprocessor");
  return params;
}

WeightGainZr::WeightGainZr(const InputParameters & parameters)
  : GeneralPostprocessor(parameters),
    _wg(0),
    _temperature(getParam<Real>("temperature")),
    _flux_integral(getPostprocessorValue("flux_integral"))
{
}

void
WeightGainZr::initialize()
{
}

void
WeightGainZr::execute()
{
  const Real Mo = 15.99;    // Molar mass of oxygen
  const Real Na(6.022e23);

  //"Initial" weight gain is temperature dependent
  Real wg0 = 2.7829;
  if (_temperature == 1273.15)
  {
    wg0 = 1.0540;
  }
  else if (_temperature == 1373.15)
  {
    wg0 = 1.7734;
  }
  else if (_temperature == 1473.15)
  {
    wg0 = 2.7829;
  }
  else if (_temperature == 1573.15)
  {
    wg0 = 4.1255;
  }
  else if (_temperature == 1673.15)
  {
    wg0 = 6.1440;
  }
  else if (_temperature == 1773.15)
  {
    wg0 = 8.2649;
  }
  else
  {
    wg0 = 6.0413e-3 * exp(4.1196e-3 * _temperature);
  }
 _wg = wg0 + 0.1 * Mo / Na * _flux_integral;  ///Weight gain in mg/cm²
}

Real
WeightGainZr::getValue()
{
    return _wg;
}
