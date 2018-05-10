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
class LineSegmentCutSetUserObject;

template <>
InputParameters validParams<LineSegmentLevelSetAux>();

/**
 * Calculate level set values for an interface that is defined by a set of line segments
 */
class LineSegmentLevelSetAux : public AuxKernel
{
public:
  LineSegmentLevelSetAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;
  virtual void compute() override;

  Real calculateSignedDistance(Point p);

  const LineSegmentCutSetUserObject * _linesegment_uo;

  std::vector<Real> _cut_data;
};

#endif // LINESEGMENTLEVELSETAUX_H
