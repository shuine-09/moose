//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "C4DiffusionCoefBeta.h"
#include "MooseUtils.h"

registerMooseObject("MooseApp", C4DiffusionCoefBeta);
registerMooseObject("MooseApp", ADC4DiffusionCoefBeta);

template <bool is_ad>
InputParameters
C4DiffusionCoefBetaTempl<is_ad>::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredParam<std::vector<std::string>>(
      "prop_names", "The names of the properties this material will have");
  //params.addRequiredParam<std::vector<Real>>("prop_values",
  //                                           "The values associated with the named properties");
  params.addParam<Real>("temperature",1473.15,"Temperature [K]. Currently only works with"
                        "homogeneous temperature in each phase.");
  params.set<MooseEnum>("constant_on") = "SUBDOMAIN";
  //params.declareControllable("prop_values");
  return params;
}

template <bool is_ad>
C4DiffusionCoefBetaTempl<is_ad>::C4DiffusionCoefBetaTempl(
    const InputParameters & parameters)
  : Material(parameters),
    _prop_names(getParam<std::vector<std::string>>("prop_names")),
    _temperature(getParam<Real>("temperature"))
    //_prop_values(getParam<std::vector<Real>>("prop_values"))
{
  unsigned int num_names = _prop_names.size();
/*
  unsigned int num_values = _prop_values.size();

  if (num_names != num_values)
    mooseError(
        "Number of prop_names must match the number of prop_values for a C4DiffusionCoefBeta!");
*/
  _num_props = num_names;

  _properties.resize(num_names);

  _prop_values.resize(num_names);

  for (unsigned int i = 0; i < _num_props; i++)
    _properties[i] = &declareGenericProperty<Real, is_ad>(_prop_names[i]);


  Real diffusivity_beta = 60 ;
  if (MooseUtils::absoluteFuzzyEqual(_temperature,633.15,1))
  {
    diffusivity_beta = 4.8e-4;
  }
  else if (MooseUtils::absoluteFuzzyEqual(_temperature,1223.15,1))
  {
    diffusivity_beta = 24;
  }
  else if (MooseUtils::absoluteFuzzyEqual(_temperature,1273.15,1))
  {
    diffusivity_beta = 38;
  }
  else if (MooseUtils::absoluteFuzzyEqual(_temperature,1373.15,1))
  {
    diffusivity_beta = 37;
  }
  else if (MooseUtils::absoluteFuzzyEqual(_temperature,1473.15,1))
  {
    diffusivity_beta = 60;
  }
  else if (MooseUtils::absoluteFuzzyEqual(_temperature,1573.15,1))
  {
    diffusivity_beta = 201.07;
  }
  else if (MooseUtils::absoluteFuzzyEqual(_temperature,1673.15,1))
  {
    diffusivity_beta = 856.19;
  }
  else if (MooseUtils::absoluteFuzzyEqual(_temperature,1773.15,1))
  {
    diffusivity_beta = 927.85;
  }
  else
  {
    diffusivity_beta = 0.110 * exp(-33391/1.987/_temperature) * 1e8;
  }

  _prop_values[0] = diffusivity_beta;
}

template <bool is_ad>
void
C4DiffusionCoefBetaTempl<is_ad>::initQpStatefulProperties()
{
  computeQpProperties();
}

template <bool is_ad>
void
C4DiffusionCoefBetaTempl<is_ad>::computeQpProperties()
{
  for (unsigned int i = 0; i < _num_props; i++)
    (*_properties[i])[_qp] = _prop_values[i];
}

template class C4DiffusionCoefBetaTempl<false>;
template class C4DiffusionCoefBetaTempl<true>;
