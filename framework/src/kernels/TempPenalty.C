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

#include "TempPenalty.h"


template<>
InputParameters validParams<TempPenalty>()
{
  InputParameters params = validParams<Kernel>();
  params.addParam<Real>("velocity", 0.0, "Plate descent speed.");
  params.addParam<Real>("diffusion", 0.0, "Thermal diffusion constant.");
  params.addParam<Real>("cold_temperature", 0.0, "Cold bath temperature.");
  params.addParam<Real>("delta_temperature", 0.0, "Temperature gap.");
  params.addParam<Real>("cold_boundary", 0.0, "Cold temperature boundary.");
  params.addParam<Real>("hot_boundary", 0.0, "Hot temperatrue boundary.");
  params.addParam<Real>("alpha", 1000, "penalty parameters.");
  return params;
}

TempPenalty::TempPenalty(const InputParameters & parameters) :
    Kernel(parameters),
    _velocity(getParam<Real>("velocity")),
    _diffusion(getParam<Real>("diffusion")),
    _cold_temp(getParam<Real>("cold_temperature")),
    _delta_temp(getParam<Real>("delta_temperature")),
    _cold_bdry(getParam<Real>("cold_boundary")),
    _hot_bdry(getParam<Real>("hot_boundary")),
    _alpha(getParam<Real>("alpha"))
{
}

TempPenalty::~TempPenalty()
{
}

Real
TempPenalty::computeQpResidual()
{

  Real cf = _t * _velocity + _cold_bdry;
  Real hf = _t * _velocity + _hot_bdry;
  Real x = _q_point[_qp](0);

  Real u0 = 0.0;
  if (x <= cf)
  {
    u0 = _cold_temp;
    return _alpha * (_u[_qp] - u0) * _test[_i][_qp];
  }
  else if (x >= hf)
  {
    u0 = _cold_temp + _delta_temp;
    return _alpha * (_u[_qp] - u0) * _test[_i][_qp];
  }
  else
  {
    return 0.0;
  }
}

Real
TempPenalty::computeQpJacobian()
{
  Real cf = _t * _velocity + _cold_bdry;
  Real hf = _t * _velocity + _hot_bdry;
  Real x = _q_point[_qp](0);

  Real u0 = 0.0;
  if (x <= cf)
  {
    u0 = _cold_temp;
    return _alpha * _phi[_j][_qp] * _test[_i][_qp];
  }
  else if (x >= hf)
  {
    u0 = _cold_temp + _delta_temp;
    return _alpha * _phi[_j][_qp] * _test[_i][_qp];
  }
  else
  {
    return 0.0;
  }
}
