//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "OxideLayerThickness.h"

registerMooseObject("XFEMApp", OxideLayerThickness);

template <>
InputParameters
validParams<OxideLayerThickness>()
{
  InputParameters params = validParams<ElementPostprocessor>();
  params.addRequiredParam<UserObjectName>("moving_line_segments",
                                          "The MovingLineSegmentCutSetUserObject user object name");
  params.addParam<unsigned int>("cut_data_index", 0, "The index of the cut data");
  return params;
}

OxideLayerThickness::OxideLayerThickness(const InputParameters & parameters)
  : ElementPostprocessor(parameters),
    _moving_line_segments(
        &getUserObject<MovingLineSegmentCutSetUserObject>("moving_line_segments")),
    _cut_data_index(getParam<unsigned int>("cut_data_index"))
{
}

Real
OxideLayerThickness::getValue()
{
  return 6e-4 - _moving_line_segments->getCutDataXCoord(_cut_data_index);
}
