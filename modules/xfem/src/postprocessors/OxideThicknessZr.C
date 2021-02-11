//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "OxideThicknessZr.h"

registerMooseObject("MooseApp", OxideThicknessZr);

defineLegacyParams(OxideThicknessZr);

InputParameters
OxideThicknessZr::validParams()
{
  InputParameters params = GeneralPostprocessor::validParams();
  params.addParam<PostprocessorName>("oxide_alpha_pos", "The name of the postprocessor giving the oxide/alpha interface position");
  return params;
}

OxideThicknessZr::OxideThicknessZr(const InputParameters & parameters)
  : GeneralPostprocessor(parameters),
    _delta(0),
    _x_ox_a(getPostprocessorValue("oxide_alpha_pos"))
{
}

void
OxideThicknessZr::initialize()
{
}

void
OxideThicknessZr::execute()
{
  const Real PBR(1.55); //Zr Pilling-Bedworth ratio
 _delta = (600 - _x_ox_a) * PBR;  //Oxide thickness [um]
}

Real
OxideThicknessZr::getValue()
{
    return _delta;
}
