//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef OXIDELAYERTHICKNESS_H
#define OXIDELAYERTHICKNESS_H

#include "ElementPostprocessor.h"
#include "MovingLineSegmentCutSetUserObject.h"
// Forward Declarations
class OxideLayerThickness;

template <>
InputParameters validParams<OxideLayerThickness>();

class OxideLayerThickness : public ElementPostprocessor
{
public:
  OxideLayerThickness(const InputParameters & parameters);

  void initialize() override {}
  void execute() override{};
  void finalize() override{};
  void threadJoin(const UserObject & user_object) override{};

  virtual Real getValue() override;

protected:
  const MovingLineSegmentCutSetUserObject * const _moving_line_segments;
  const unsigned int _cut_data_index;
};

#endif // OXIDELAYERTHICKNESS_H
