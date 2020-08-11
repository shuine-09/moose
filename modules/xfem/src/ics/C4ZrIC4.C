//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "C4ZrIC4.h"
#include "MooseUtils.h"

registerMooseObject("MooseApp", C4ZrIC4);

defineLegacyParams(C4ZrIC4);

InputParameters
C4ZrIC4::validParams()
{
  InputParameters params = InitialCondition::validParams();
  //params.addRequiredParam<FunctionName>("function", "The initial condition function.");
  params.addRequiredParam<Real>("temperature", "The temperature of the cladding (homogeneous temperature only)");
  params.addClassDescription("An initial condition for the C4 model for high temperature "
                            "corrosion of Zircaloy-4");
  return params;
}

C4ZrIC4::C4ZrIC4(const InputParameters & parameters)
  : InitialCondition(parameters),_temperature(getParam<Real>("temperature"))
{
  if (MooseUtils::absoluteFuzzyEqual(_temperature,1273.15,1))
  {
    _x_a_b = 591.4;
    _x_ox_a = 595.6;
    _x_b_break = 549.645326;
    _x_a_break = 593.9469;
    _C_a_break = 0.0571;
  }
  else if (MooseUtils::absoluteFuzzyEqual(_temperature,1373.15,1))
  {
    _x_a_b = 583.6;
    _x_ox_a = 593.2;
    _x_b_break = 547.657598;
    _x_a_break = 588.6630;
    _C_a_break = 0.0772;
  }
  else if (MooseUtils::absoluteFuzzyEqual(_temperature,1473.15,1))
  {
    _x_a_b = 572.2;
    _x_ox_a = 590.0;
    _x_b_break = 530.82328;
    _x_a_break = 580.7314;
    _C_a_break = 0.0907;
  }
  else if (MooseUtils::absoluteFuzzyEqual(_temperature,1573.15,1))
  {
    _x_a_b = 560.7;
    _x_ox_a = 586.0;
    _x_b_break = 481.08946;
    _x_a_break = 572.0145;
    _C_a_break = 0.1026;
  }
  else if (MooseUtils::absoluteFuzzyEqual(_temperature,1673.15,1))
  {
    _x_a_b = 539.5;
    _x_ox_a = 582.0;
    _x_b_break = 370.895773;
    _x_a_break = 556.2216;
    _C_a_break = 0.1141;
  }
  else if (MooseUtils::absoluteFuzzyEqual(_temperature,1773.15,1))
  {
    _x_a_b = 516.3;
    _x_ox_a = 576.2;
    _x_b_break = 358.985546;
    _x_a_break = 539.5872;
    _C_a_break = 0.1258;
  }
  else
  {
    _x_a_b = -1.9259*1e-7*pow(_temperature,3) + 6.7254*1e-4*pow(_temperature,2) - 0.84697*_temperature + 977.01  ;
    _x_ox_a = -3.6071*1e-5*pow(_temperature,2) + 7.1427*1e-2*_temperature + 563.11  ;
    _x_b_break = 3.6833*1e-8*pow(_temperature,4) - 2.1993*1e-4*pow(_temperature,3) + 0.48902*pow(_temperature,2) - 480.24*_temperature + 176420  ;
    _x_a_break = -1.5750*1e-4*pow(_temperature,2) + 0.37183*_temperature + 375.50;
    _C_a_break = 3.3876*1e-10*pow(_temperature,3) - 1.6378*1e-6*pow(_temperature,2) + 2.7475*1e-3*_temperature - 1.4851;
  }

  const Real Zr_PBR = 1.55;
  const Real xo_ox = 0.66666666667;
  _C_ox = 3/Zr_PBR * xo_ox;

  const Real xo_ox_a = 0.667118123 - 1.10606e-5 * _temperature;
  _C_ox_a = 3/Zr_PBR * xo_ox_a;

  Real xo_a_ox;
  if (_temperature > 473.15 && _temperature < 1478.15)
  {
    xo_a_ox = (28.6 + exp(-6748/_temperature + 4.748)) * 1e-2;
  }
  else if (_temperature > 1478.15 && _temperature < 1798.15)
  {
    xo_a_ox = (28.6 + exp(-6301/_temperature + 4.460)) * 1e-2;
  }
  else if (_temperature > 1798.15 && _temperature < 2338.15)
  {
    xo_a_ox = (28.6 + exp(-7012/_temperature + 8.434 - 3.521e-3 * _temperature )) * 1e-2;
  }
  else
  {
    xo_a_ox = 28.6 *1e-2;
  }
  _C_a_ox = xo_a_ox/(1-xo_a_ox);

  Real xo_b_a = (9.59e-3 * (_temperature - 1136) + 4.72e-6 * pow(_temperature - 1136,2) - 4.35e-9 * pow(_temperature - 1136,3)) * 1e-2;
  Real xo_a_b = (45.86e-3 * (_temperature - 1136) - 44.77e-6 * pow(_temperature - 1136,2) + 17.40e-9 * pow(_temperature - 1136,3)) * 1e-2; //the original one, not the weak equivalent
  _C_a_b = xo_a_b / (1 - xo_a_b);
  _C_b_a = xo_b_a / (1 - xo_b_a);

  _C_ox_a_weak = _C_ox_a - (_C_a_b - _C_b_a);
  _C_ox_weak = _C_ox - (_C_ox_a - _C_a_ox) - (_C_a_b - _C_b_a);

  _grad_ba = (_C_b_a - _C_b) / (_x_a_b - _x_b_break) ;
  _grad_ab = (_C_a_break - _C_b_a) / (_x_a_break - _x_a_b);
  _grad_aox = (_C_ox_a_weak - _C_a_break) / (_x_ox_a - _x_a_break);
  _grad_ox = (_C_ox_weak - _C_ox_a_weak) / (600 - _x_ox_a);

}

Real
C4ZrIC4::value(const Point & p)
{
  if (p(0)<_x_b_break)
  {
    return _C_b;
  }
  else if (p(0)<_x_a_b)
  {
    return (_C_b_a + (p(0)-_x_a_b)*_grad_ba);
  }
  else if (p(0)<_x_a_break)
  {
    return (_C_b_a + (p(0)-_x_a_b)*_grad_ab);
  }
  else if (p(0)<_x_ox_a)
  {
    return (_C_ox_a_weak + (p(0)-_x_ox_a)*_grad_aox);
  }
  else
  {
    return (_C_ox_weak + (p(0)-600) * _grad_ox);
  }
}

/**RealGradient
C4ZrIC4::gradient(const Point & p)
{
  return _func.gradient(_t, p);
}
*/
