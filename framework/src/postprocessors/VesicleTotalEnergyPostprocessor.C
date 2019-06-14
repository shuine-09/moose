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

#include "VesicleTotalEnergyPostprocessor.h"

registerMooseObject("MooseApp", VesicleTotalEnergyPostprocessor);

template <>
InputParameters
validParams<VesicleTotalEnergyPostprocessor>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addRequiredCoupledVar("variable",
                               "The name of the variable that this userobject applies to");
  params.addParam<Real>("spontaneous_curvature", 0.0, "Spontaneous of vesicle.");
  params.addParam<Real>("epsilon", 0.01, "The interfacial penalty parameter.");
  params.addParam<bool>("rz", false, "The RZ coordindate system.");
  return params;
}

VesicleTotalEnergyPostprocessor::VesicleTotalEnergyPostprocessor(const InputParameters & parameters)
  : ElementIntegralPostprocessor(parameters),
    MooseVariableInterface<Real>(this, false),
    _u(coupledValue("variable")),
    _grad_u(coupledGradient("variable")),
    _second_u(coupledSecond("variable")),
    _C(getParam<Real>("spontaneous_curvature")),
    _epsilon(getParam<Real>("epsilon")),
    _rz(getParam<bool>("rz")),
    _qp(0)
{
}

void
VesicleTotalEnergyPostprocessor::threadJoin(const UserObject & y)
{
  const VesicleTotalEnergyPostprocessor & pps =
      static_cast<const VesicleTotalEnergyPostprocessor &>(y);
  _integral_value += pps._integral_value;
}

Real
VesicleTotalEnergyPostprocessor::computeQpIntegral()
{

  Real rz_coord = _q_point[_qp](0);

  Real lap_u = _second_u[_qp].tr();

  if (_rz)
    lap_u += _grad_u[_qp](0) / rz_coord;

  //  Real value =
  //      0.5 * pow(lap_u - 1.0 / (_epsilon * _epsilon) * (_u[_qp] * _u[_qp] - 1) * _u[_qp] -
  //      _C, 2.0);
  Real value = 0.5 * pow(lap_u - 1.0 / (_epsilon * _epsilon) * (_u[_qp] * _u[_qp] - 1) *
                                     (_u[_qp] + _C * _epsilon),
                         2.0);

  return value;
}
