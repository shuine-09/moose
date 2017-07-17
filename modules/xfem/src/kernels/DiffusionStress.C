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

#include "DiffusionStress.h"

template <>
InputParameters
validParams<DiffusionStress>()
{
  InputParameters params = validParams<Kernel>();
  params.addCoupledVar("pressure", "The pressure aux variable name");
  params.addParam<Real>("coefficent", 1, "Coefficient of pressure term");
  params.addClassDescription("The Laplacian operator ($-\\nabla \\cdot \\nabla u$), with the weak "
                             "form of $(\\nabla \\phi_i, \\nabla u_h)$.");
  return params;
}

DiffusionStress::DiffusionStress(const InputParameters & parameters)
  : Kernel(parameters),
    _pressure(coupledValue("pressure")),
    _grad_pressure(coupledGradient("pressure")),
    _coef(getParam<Real>("coefficent"))
{
}

Real
DiffusionStress::computeQpResidual()
{
  return _grad_u[_qp] * _grad_test[_i][_qp] -
         _coef * _u[_qp] * _grad_test[_i][_qp] * _grad_pressure[_qp];
}

Real
DiffusionStress::computeQpJacobian()
{
  return _grad_phi[_j][_qp] * _grad_test[_i][_qp] -
         _coef * _phi[_j][_qp] * _grad_test[_i][_qp] * _grad_pressure[_qp];
}
