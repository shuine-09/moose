/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "VesicleIC.h"

registerMooseObject("PhaseFieldApp", VesicleIC);

template <>
InputParameters
validParams<VesicleIC>()
{
  InputParameters params = validParams<InitialCondition>();
  params.addParam<Point>("center", Point(0, 0, 0), "The center coordinate");
  params.addParam<Real>("major", "The major axis of the ellipse.");
  params.addParam<Real>("minor", "The minor axis of the ellipse.");
  params.addParam<Real>("epsilon", 0.01, "The interfacial penalty parameter.");
  return params;
}

VesicleIC::VesicleIC(const InputParameters & parameters)
  : InitialCondition(parameters),
    _center(parameters.get<Point>("center")),
    _major_axis(getParam<Real>("major")),
    _minor_axis(getParam<Real>("minor")),
    _epsilon(getParam<Real>("epsilon"))
{
}

Real
VesicleIC::value(const Point & p)
{
  Real value = 0.0;

  value =
      pow(pow((p(0) - _center(0)) / _major_axis, 2.0) + pow((p(1) - _center(1)) / _minor_axis, 2.0),
          0.5) -
      1.0;

  value = std::tanh(value / sqrt(2.0) / _epsilon);

  return value;
}
