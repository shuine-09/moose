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
#include "WeibullMaterial.h"

registerMooseObject("TensorMechanicsApp", WeibullMaterial);

template <>
InputParameters
validParams<WeibullMaterial>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<Real>(
      "weibull_modulus", "The Weibull modulus quantifying observed variability in the property.");
  params.addRequiredParam<Real>(
      "specimen_material_property",
      "The median value of material property observed in the laboratory.");
  params.addRequiredParam<Real>("specimen_volume", "Specimen volume used in the laboratory.");
  params.addParam<Real>("high", 0.0, "Specimen volume used in the laboratory.");
  params.addParam<Real>("low", 0.0,"Specimen volume used in the laboratory.");
  return params;
}

WeibullMaterial::WeibullMaterial(const InputParameters & parameters)
  : Material(parameters),
    _weibull(declareProperty<Real>("gc_prop")),
    _weibull_old(getMaterialPropertyOld<Real>("gc_prop")),
    _weibull_modulus(getParam<Real>("weibull_modulus")),
    _specimen_volume(getParam<Real>("specimen_volume")),
    _specimen_material_property(getParam<Real>("specimen_material_property")),
    _high(getParam<Real>("high")),
    _low(getParam<Real>("low"))
{
  // Setup the random number generation
  setRandomResetFrequency(EXEC_INITIAL);
}

void
WeibullMaterial::initQpStatefulProperties()
{
  if (_qp == 0)
  {
    if (std::abs(_weibull_modulus) > 1.0e-5)
    {
      Real rn = getRandomReal();
      if (_high > 1.0)
      {
        Point p = _q_point[_qp];
        Real x = p(0);
        Real z = p(2);
        Real r = std::sqrt(x*x+z*z);
        Real mp = (_high - _low) / 0.0041*r + _low;
      _eta = mp *
             std::pow(_specimen_volume * std::log(rn) / (_current_elem->volume() * std::log(0.5)),
                      1.0 / (Real)_weibull_modulus);
      }
      else
      {
      _eta = _specimen_material_property *
             std::pow(_specimen_volume * std::log(rn) / (_current_elem->volume() * std::log(0.5)),
                      1.0 / (Real)_weibull_modulus);
      }
    }
    _weibull[_qp] = _eta;
  }
  else
  {
    _weibull[_qp] = _eta;
  }
}

void
WeibullMaterial::computeQpProperties()
{
}
