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

#include "HeatTestAux.h"
#include "Function.h"

template<>
InputParameters validParams<HeatTestAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addParam<Real>("velocity", 0.0, "Plate descent speed.");
  params.addParam<Real>("diffusion", 0.0, "Thermal diffusion constant.");
  params.addParam<Real>("cold_temperature", 0.0, "Cold bath temperature.");
  params.addParam<Real>("delta_temperature", 0.0, "Temperature gap.");
  params.addParam<Real>("cold_boundary", 0.0, "Cold temperature boundary.");
  params.addParam<Real>("hot_boundary", 0.0, "Hot temperatrue boundary.");
  return params;
}

HeatTestAux::HeatTestAux(const InputParameters & parameters) :
    AuxKernel(parameters),
    _velocity(getParam<Real>("velocity")),
    _diffusion(getParam<Real>("diffusion")),
    _cold_temp(getParam<Real>("cold_temperature")),
    _delta_temp(getParam<Real>("delta_temperature")),
    _cold_bdry(getParam<Real>("cold_boundary")),
    _hot_bdry(getParam<Real>("hot_boundary"))
{
}

Real
HeatTestAux::computeValue()
{
  Real l = -1.0 * (_t * _velocity + _cold_bdry);
  Real P = _diffusion / _velocity;
  Real x = (*_current_node)(0);

  if (isNodal())
  {
    return _cold_temp + _delta_temp * (1 - std::exp(-(x+l)/P)) * step_func(x+l);
  }
  else
    mooseError("HeatTestAux only supports Nodal AuxVariable");

  /*
  Real cb = (_t * _velocity + _cold_bdry);
  Real hb = (_t * _velocity + _hot_bdry);
  Real x = (*_current_node)(0);

  if (isNodal())
  {
    if (x < cb)
      return _cold_temp;
    else if (x >= cb && x <= hb)
      return (x-cb)/(hb-cb)*(_delta_temp) + _cold_temp;
    else
      return _delta_temp + _cold_temp;
  }
  else
    mooseError("HeatTestAux only supports Nodal AuxVariable");
  */

}

Real
HeatTestAux::step_func(Real x)
{
  if (std::fabs(x) < 1.0e-10)
    return 0.5;
  else if (x > 0.0)
    return 1.0;
  else return 0.0;
}
