//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MeltPoolHeight.h"

registerMooseObject("LevelSetApp", MeltPoolHeight);

defineLegacyParams(MeltPoolHeight);

InputParameters
MeltPoolHeight::validParams()
{
  InputParameters params = ElementPostprocessor::validParams();
  params.addClassDescription("Compute melt pool height.");
  params.addRequiredCoupledVar("level_set", "Level set variable");
  return params;
}

MeltPoolHeight::MeltPoolHeight(const InputParameters & parameters)
  : ElementPostprocessor(parameters),
    _ls(coupledValue("level_set")),
    _f_l(getMaterialProperty<Real>("liquid_mass_fraction")),
    _min_value(std::numeric_limits<Real>::max()),
    _max_value(-std::numeric_limits<Real>::max())
{
}

void
MeltPoolHeight::initialize()
{
  _max_value = -std::numeric_limits<Real>::max(); // start w/ the min
  _min_value = std::numeric_limits<Real>::max();  // start w/ the max
}

void
MeltPoolHeight::execute()
{
  for (_qp = 0; _qp < _qrule->n_points(); _qp++)
    computeQpValue();
}

void
MeltPoolHeight::computeQpValue()
{
  if (_ls[_qp] < 0.5 && _f_l[_qp] > 0.99)
  {
    _max_value = std::max(_max_value, _q_point[_qp](1));

    _min_value = std::min(_min_value, _q_point[_qp](1));
  }
}

void
MeltPoolHeight::threadJoin(const UserObject & uo)
{
  const MeltPoolHeight & pps = static_cast<const MeltPoolHeight &>(uo);

  _max_value = std::max(_max_value, pps._max_value);

  _min_value = std::min(_min_value, pps._min_value);
}

Real
MeltPoolHeight::getValue()
{
  gatherMax(_max_value);

  gatherMin(_min_value);

  _height = _max_value - _min_value;

  return _height < 0 ? 0 : _height;
}
