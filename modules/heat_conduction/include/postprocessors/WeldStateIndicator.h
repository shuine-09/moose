//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef WELDSTATEINDICATOR_H
#define WELDSTATEINDICATOR_H

#include "ElementUserObject.h"

// Forward Declarations
class WeldStateIndicator;

// Input parameters
template <>
InputParameters validParams<WeldStateIndicator>();

/// A postprocessor for collecting the nodal min or max value
class WeldStateIndicator : public ElementUserObject
{
public:
  enum WeldStateType
  {
    BEFORE,
    HEATING,
    COOLING,
    AFTER
  };

  /**
   * Class constructor
   * @param parameters The input parameters
   */
  WeldStateIndicator(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void threadJoin(const UserObject & y) override;
  virtual void execute() override;
  virtual void finalize() override;

  WeldStateType getWeldState() const { return _state; };

protected:
  /// Holds the solution at current quadrature points
  const VariableValue & _u;

  Real _min_u;

  bool _has_melted;

  WeldStateType _state;

  Real _start_time;
  Real _end_time;
  Real _melting_temp;
};

#endif
