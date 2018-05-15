//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef MOVINGMovingLineSegmentCutSetUserObject_H
#define MOVINGMovingLineSegmentCutSetUserObject_H

#include "GeometricCut2DUserObject.h"
#include "ElementPairLocator.h"
#include "GeneralVectorPostprocessor.h"
#include "VectorPostprocessorInterface.h"

// Forward declarations
class MovingLineSegmentCutSetUserObject;
class PointValueAtXFEMInterface;

template <>
InputParameters validParams<MovingLineSegmentCutSetUserObject>();

class MovingLineSegmentCutSetUserObject : public GeometricCut2DUserObject,
                                          protected VectorPostprocessorInterface

{
public:
  MovingLineSegmentCutSetUserObject(const InputParameters & parameters);

  virtual void initialize() override;

  virtual void execute() override;

  virtual void finalize() override;

  virtual const std::vector<Point>
  getCrackFrontPoints(unsigned int num_crack_front_points) const override;

  Real getLocationX() const;

  virtual std::vector<Real> getCutData() const { return _cut_data; };

protected:
  const PointValueAtXFEMInterface * _interface_value_uo;

  std::vector<Real> _cut_data;

  /// The variable number of the solution variable we using to calcuate velocity
  const unsigned int _var_number;

  /// system reference
  const System & _system;

  /// the subproblem solution vector
  const NumericVector<Number> * _solution;
};

#endif // MOVINGMovingLineSegmentCutSetUserObject_H
