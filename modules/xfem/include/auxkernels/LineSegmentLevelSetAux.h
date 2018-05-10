//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef LINESEGMENTLEVELSETAUX_H
#define LINESEGMENTLEVELSETAUX_H

#include "AuxKernel.h"

// Forward Declarations
class LineSegmentLevelSetAux;
class MovingLineSegmentCutSetUserObject;

template <>
InputParameters validParams<LineSegmentLevelSetAux>();

/**
 * Function auxiliary value
 */
class LineSegmentLevelSetAux : public AuxKernel
{
public:
  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  LineSegmentLevelSetAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;
  virtual void compute() override;

  Real calculateSignedDistance(Point p);

  /// Function being used to compute the value of this kernel
  const MovingLineSegmentCutSetUserObject * _linesegment_uo;

  std::vector<Real> _cut_data;
};

#endif // LINESEGMENTLEVELSETAUX_H
