/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "XFEMMaterialTensorMarkerUserObject.h"
#include "RankTwoScalarTools.h"
#include "libmesh/quadrature.h"

template <>
InputParameters
validParams<XFEMMaterialTensorMarkerUserObject>()
{
  InputParameters params = validParams<XFEMMarkerUserObject>();
  params.addParam<MooseEnum>(
      "scalar_type", RankTwoScalarTools::scalarOptions(), "Type of scalar output");
  params.addRequiredParam<MaterialPropertyName>("rank_two_tensor",
                                                "The rank two material tensor name");
  params.addRequiredParam<Real>("threshold", "The threshold for crack growth.");
  params.addRequiredParam<bool>(
      "average", "Should the tensor quantity be averaged over the quadruature points?");
  params.addParam<Real>(
      "random_range", 0.0, "Range of a uniform random distribution for the threshold");
  params.addParam<Point>(
      "point1",
      Point(0, 0, 0),
      "Start point for axis used to calculate some cylinderical material tensor quantities");
  params.addParam<Point>("point2",
                         Point(0, 1, 0),
                         "End point for axis used to calculate some material tensor quantities");
  params.addParam<bool>("use_weibull", false, "Use weibull distribution for material strength");
  return params;
}

XFEMMaterialTensorMarkerUserObject::XFEMMaterialTensorMarkerUserObject(
    const InputParameters & parameters)
  : XFEMMarkerUserObject(parameters),
    _use_weibull(getParam<bool>("use_weibull")),
    _tensor(getMaterialProperty<RankTwoTensor>("rank_two_tensor")),
    _scalar_type(getParam<MooseEnum>("scalar_type")),
    _point1(parameters.get<Point>("point1")),
    _point2(parameters.get<Point>("point2")),
    _threshold(getParam<Real>("threshold")),
    _average(getParam<bool>("average")),
    _random_range(getParam<Real>("random_range")),
    _weibull(getMaterialProperty<Real>("weibull"))
{
  setRandomResetFrequency(EXEC_INITIAL);
}

bool
XFEMMaterialTensorMarkerUserObject::doesElementCrack(Point & direction)
{
  bool does_it_crack = false;
  unsigned int numqp = _qrule->n_points();

  Real rnd_mult = (1.0 - _random_range / 2.0) + _random_range * getRandomReal();

  Real perturbed_threshold = 0.0;
  if (_use_weibull) // use weibull
    perturbed_threshold = _threshold * _weibull[0];
  else
    perturbed_threshold = _threshold * rnd_mult;

  Real average_quantity = 0;
  if (_average)
  {
    RankTwoTensor average_tensor;
    average_tensor.zero();
    for (unsigned int qp = 0; qp < numqp; ++qp)
    {
      average_tensor += _tensor[qp];
      average_quantity += RankTwoScalarTools::getQuantity(
          _tensor[qp], _scalar_type, _point1, _point2, _q_point[qp], direction);
    }

    average_tensor *= 1.0 / (Real)numqp;
    average_quantity *= 1.0 / (Real)numqp;

    if (average_quantity > perturbed_threshold)
      does_it_crack = true;
  }
  else
  {
    unsigned int max_index = 999999;
    std::vector<Real> tensor_quantities;
    tensor_quantities.reserve(numqp);
    Real max_quantity = 0;
    std::vector<Point> directions;
    directions.resize(numqp);
    for (unsigned int qp = 0; qp < numqp; ++qp)
    {
      tensor_quantities[qp] = RankTwoScalarTools::getQuantity(
          _tensor[qp], _scalar_type, _point1, _point2, _q_point[qp], directions[qp]);
      if (directions[qp](0) == 0 && directions[qp](1) == 0 && directions[qp](2) == 0)
      {
        mooseError("Direction has zero length in XFEMMaterialTensorMarkerUserObject");
      }
      if (tensor_quantities[qp] > max_quantity)
      {
        max_quantity = tensor_quantities[qp];
        max_index = qp;
      }
    }
    if (max_quantity > perturbed_threshold)
    {
      does_it_crack = true;
      direction = directions[max_index];
    }
  }

  return does_it_crack;
}
