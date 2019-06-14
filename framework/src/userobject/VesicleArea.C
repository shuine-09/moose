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

#include "VesicleArea.h"

// libmesh includes
#include "libmesh/quadrature.h"

registerMooseObject("MooseApp", VesicleArea);

template <>
InputParameters
validParams<VesicleArea>()
{
  InputParameters params = validParams<ElementUserObject>();
  params.addRequiredCoupledVar("variable",
                               "The name of the variable that this userobject applies to");
  params.addParam<Real>("epsilon", 0.01, "The interfacial penalty parameter.");
  return params;
}

VesicleArea::VesicleArea(const InputParameters & parameters)
  : ElementUserObject(parameters),
    MooseVariableInterface<Real>(this, false),
    _u(coupledValue("variable")),
    _grad_u(coupledGradient("variable")),
    _epsilon(getParam<Real>("epsilon")),
    _qp(0),
    _integral_value(0)
{
}

void
VesicleArea::initialize()
{
  _integral_value = 0;
}

void
VesicleArea::execute()
{
  _integral_value += computeIntegral();
}

void
VesicleArea::finalize()
{
  gatherSum(_integral_value);
}

Real
VesicleArea::getValue() const
{
  return _integral_value;
}

void
VesicleArea::threadJoin(const UserObject & y)
{
  const VesicleArea & pps = static_cast<const VesicleArea &>(y);
  _integral_value += pps._integral_value;
}

Real
VesicleArea::computeIntegral()
{
  Real sum = 0;

  for (_qp = 0; _qp < _qrule->n_points(); _qp++)
    sum += _JxW[_qp] * _coord[_qp] * computeQpIntegral();
  return sum;
}

Real
VesicleArea::computeQpIntegral()
{
  Real value = 0.5 * _epsilon * _grad_u[_qp] * _grad_u[_qp] +
               1.0 / (4.0 * _epsilon) * pow((_u[_qp] * _u[_qp] - 1.0), 2.0);
  value *= 3.0 / 2.0 / std::sqrt(2.0);
  return value;
}
