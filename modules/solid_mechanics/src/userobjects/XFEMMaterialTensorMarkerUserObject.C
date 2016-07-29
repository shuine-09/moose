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

#include "XFEMMaterialTensorMarkerUserObject.h"

template<>
InputParameters validParams<XFEMMaterialTensorMarkerUserObject>()
{
  InputParameters params = validParams<XFEMMarkerUserObject>();
  params += validParams<MaterialTensorCalculator>();
  params.addRequiredParam<std::string>("tensor", "The material tensor name.");
  params.addRequiredParam<Real>("threshold", "The threshold for crack growth.");
  params.addRequiredParam<bool>("average", "Should the tensor quantity be averaged over the quadruature points?");
  params.addParam<Real>("random_range",0.0,"Range of a uniform random distribution for the threshold");
  params.addParam<bool>("use_weibull", false,"Use weibull distribution to initiate crack?");
  //params.addCoupledVar("weibull","Weibull aux variable");
  return params;
}

XFEMMaterialTensorMarkerUserObject::XFEMMaterialTensorMarkerUserObject(const InputParameters & parameters):
  XFEMMarkerUserObject(parameters),
  _use_weibull(getParam<bool>("use_weibull")),
  //_weibull(_use_weibull? coupledValue("weibull") : _zero),
  _material_tensor_calculator(parameters),
  _tensor(getMaterialProperty<SymmTensor>(getParam<std::string>("tensor"))),
  _threshold(getParam<Real>("threshold")),
  _average(getParam<bool>("average")),
  _random_range(getParam<Real>("random_range")),
  _weibull_eta(getMaterialProperty<Real>("weibull_eta"))
{
  setRandomResetFrequency(EXEC_INITIAL);
}

bool
XFEMMaterialTensorMarkerUserObject::doesElementCrack(RealVectorValue &direction)
{
  bool does_it_crack = false;
  unsigned int numqp = _qrule->n_points();

  Real rnd_mult = (1.0 - _random_range/2.0) + _random_range*getRandomReal();

  Real perturbed_threshold = 0.0;
  if (_use_weibull) //use weibull 
    perturbed_threshold = _weibull_eta[0];
  else
    perturbed_threshold = _threshold * rnd_mult;

  if (_average)
  {
    SymmTensor average_tensor;
    for ( unsigned int qp = 0; qp < numqp; ++qp )
    {
      average_tensor += _tensor[qp];
    }
    average_tensor *= 1.0/(Real)numqp;
    Real tensor_quantity = _material_tensor_calculator.getTensorQuantity(average_tensor,&_q_point[0],direction);
    if(tensor_quantity > perturbed_threshold)
      does_it_crack = true;
  }
  else
  {
    unsigned int max_index = 999999;
    std::vector<Real> tensor_quantities;
    tensor_quantities.reserve(numqp);
    Real max_quantity = 0;
    std::vector<RealVectorValue> directions;
    directions.resize(numqp);
    for ( unsigned int qp = 0; qp < numqp; ++qp )
    {
      tensor_quantities[qp] = _material_tensor_calculator.getTensorQuantity(_tensor[qp],&_q_point[qp],directions[qp]);
      if (directions[qp](0) == 0 &&
          directions[qp](1) == 0 &&
          directions[qp](2) == 0)
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
