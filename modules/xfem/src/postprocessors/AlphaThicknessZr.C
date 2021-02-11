//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AlphaThicknessZr.h"

registerMooseObject("MooseApp", AlphaThicknessZr);

defineLegacyParams(AlphaThicknessZr);

InputParameters
AlphaThicknessZr::validParams()
{
  InputParameters params = GeneralPostprocessor::validParams();
  params.addParam<PostprocessorName>("oxide_alpha_pos", "The name of the postprocessor giving the oxide/alpha interface position");
  params.addParam<PostprocessorName>("alpha_beta_pos", "The name of the postprocessor giving the alpha/beta interface position");
  return params;
}

AlphaThicknessZr::AlphaThicknessZr(const InputParameters & parameters)
  : GeneralPostprocessor(parameters),
    _d_alpha(0),
    _x_ox_a(getPostprocessorValue("oxide_alpha_pos")),
    _x_a_b(getPostprocessorValue("alpha_beta_pos"))
{
}

void
AlphaThicknessZr::initialize()
{
}

void
AlphaThicknessZr::execute()
{
 _d_alpha = _x_ox_a - _x_a_b; //Alpha layer thickness [um]
}

Real
AlphaThicknessZr::getValue()
{
    return _d_alpha;
}
