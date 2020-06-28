//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ElementUserObject.h"

class MeltPoolLevelSetLocation;
class Function;

/**
 *
 */
template <>
InputParameters validParams<MeltPoolLevelSetLocation>();

class MeltPoolLevelSetLocation : public ElementUserObject
{
public:
  static InputParameters validParams();

  MeltPoolLevelSetLocation(const InputParameters & parameters);

  void initialize() override;
  void execute() override;
  void threadJoin(const UserObject & uo) override;
  void finalize() override;
  const Point getLaserSpotLocation() const { return Point(_location_x, _location_y, _location_z); };

protected:
  /// Function of laser location in x coordinate
  const Function & _laser_center_x;

  /// Function of laser location in x coordinate
  const Function & _laser_center_y;

  /// Function of laser location in x coordinate
  const Function & _laser_center_z;

  /// The variable number of the level set variable we are operating on
  const unsigned int _level_set_var_number;

  /// system reference
  const System & _system;

  /// the subproblem solution vector
  const NumericVector<Number> * _solution;

  Real _location_x;
  Real _location_y;
  Real _location_z;

  Point _starting_point;
  Real _velocity;
};
