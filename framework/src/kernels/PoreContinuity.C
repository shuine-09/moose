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

#include "PoreContinuity.h"


template<>
InputParameters validParams<PoreContinuity>()
{
  InputParameters params = validParams<Kernel>();
  params.addCoupledVar("temperature", "Temperature.");
  params.addParam<bool>("debug", 0, "Whether or not to debug");
  params.addParam<Real>("alpha", 1, "alpha");
  params.addParam<Real>("beta", 1, "beta");
  params.addRequiredParam<Real>("nu", "porosity diffusion coefficient");
  return params;
}

PoreContinuity::PoreContinuity(const InputParameters & parameters) :
  Kernel(parameters),
  _has_temperature(isCoupled("temperature")),
  _T(_has_temperature ? coupledValue("temperature") : _zero),
  _grad_T(coupledGradient("temperature")),
  _debug(getParam<bool>("debug")),
  //   _pore_velocity(getMaterialProperty<RealVectorValue>("pore_velocity"))
  _pore_velocity(getMaterialProperty<Real>("pore_velocity")),
  _v_var(coupled("temperature")),
  _alpha(getParam<Real>("alpha")),
  _beta(getParam<Real>("beta")),
  _nu(getParam<Real>("nu"))
{}

PoreContinuity::~PoreContinuity()
{
}

Real
PoreContinuity::convection_scalar()
{
  // _pore_velocity is really small, like 1.e-25 initially, so maybe
  // we'll take it out of the equation (replace it with a constant value) for the moment...
  // return _alpha * _pore_velocity[_qp] / (_beta+_T[_qp]) / (_beta+_T[_qp]);
  //  return _alpha * _pore_velocity[_qp] / (_beta + _T[_qp]) / (_beta + _T[_qp]);
  return _alpha * _pore_velocity[_qp];
}

RealVectorValue
PoreContinuity::convection_velocity()
{
  // Compute the convective velocity a = (alpha * vp)/(beta + T)^2 * grad(T)
  return convection_scalar() * _grad_T[_qp];
}



Real
PoreContinuity::computeQpResidual()
{
  if (_debug)
  {
    // Moose::out << "res = " << -_u[_qp] * (convection_velocity() * _grad_test[_i][_qp]) << std::endl;
    // Moose::out << "u = " <<  _u[_qp] << std::endl;
    // Moose::out << "grad_u = " <<  _grad_u[_qp] << std::endl;
    // Moose::out << "temperature = " <<  _T[_qp] << std::endl;
    // Moose::out << "grad_T = " <<  _grad_T[_qp] << std::endl;
    // Moose::out << "grad_test = " <<  _grad_test[_i][_qp] << std::endl;
    Moose::out << "a/(b+T)^2 = " << _alpha / (_beta+_T[_qp]) / (_beta+_T[_qp]) << std::endl;
    Moose::out << "pore_velocity = " <<  _pore_velocity[_qp] << std::endl;
    Moose::out << "|a|=" << convection_velocity().norm() << std::endl;
  }

  //  return _scale * _u[_qp] / ( (1 + _T[_qp])*(1 + _T[_qp])) * _grad_T[_qp] * _grad_test[_i][_qp];// from chemotaxis example
  //  return _scale * _pore_velocity[_qp] * _u[_qp] * _grad_T[_qp] * _grad_test[_i][_qp];
  //  if (_u[_qp] < 0.5)
  //    return 0.;
  //  else

  // Compute the residual contribution -ap * grad(test) + nu * grad(p) * grad(test)
  return -_u[_qp] * (convection_velocity() * _grad_test[_i][_qp]) + _nu * (_grad_u[_qp] * _grad_test[_i][_qp]);
}

Real
PoreContinuity::computeQpJacobian()
{
  //  return _grad_phi[_j][_qp] * _grad_test[_i][_qp];// not sure of form
  //  return _scale * _pore_velocity[_qp] * _grad_test[_i][_qp];// d/du of res?
  //  if (_u[_qp] < 0.5)
  //    return 0.;
  //  else

  // Compute the jacobian contribution for -ap * grad(test) + nu * grad(p) * grad(test)
  return -_phi[_j][_qp] * (convection_velocity() * _grad_test[_i][_qp]) + _nu * (_grad_phi[_j][_qp] * _grad_test[_i][_qp]);
}

Real
PoreContinuity::computeQpOffDiagJacobian(unsigned int jvar)
{
  // Derivative with respect to temperature.  Note, we are ignoring
  // the off-diagonal contribution due to d(pore_velocity)/dT.
  if (jvar == _v_var)
    return -_u[_qp] * convection_scalar() * (_grad_phi[_j][_qp] - 2 * _phi[_j][_qp] / (_beta + _T[_qp]) * _grad_T[_qp]) * _grad_test[_i][_qp];

  return 0.0;
}
