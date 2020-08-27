//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeometricCut2DUserObject.h"

// Forward declarations

/**
*
* Added booleans to specify if the object is representing an interface in the
* weak discontinuity equivalent of the C4 model for the high-temperature
* corrosion of Zircaloy-4 (1000C to 1500C), and if so which interfaces and what
* the temperature is.
* This is used to specify the "initial" (time = 20s) positions of the _two_interfaces
* as retrieved from the finite difference Matlab implementation of the model.
*
* By defaults, all the booleans are set to false so that the class can keep
* being used to model interfaces in other applications without any change. 
*/

class LineSegmentCutSetUserObject : public GeometricCut2DUserObject
{
public:
  static InputParameters validParams();

  LineSegmentCutSetUserObject(const InputParameters & parameters);

  virtual const std::vector<Point>

  getCrackFrontPoints(unsigned int num_crack_front_points) const override;

  /**
   * Get the cut location information
   */
  virtual std::vector<Real> getCutData() const { return _cut_data; };

protected:

  std::vector<Real> _cut_data;

  bool _is_C4;

  bool _ab_interface;

  bool _oxa_interface;

  Real _temperature;
};
