//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ElementPostprocessor.h"

class MeltPoolHeight;
class Function;

/**
 *
 */
template <>
InputParameters validParams<MeltPoolHeight>();

class MeltPoolHeight : public ElementPostprocessor
{
public:
  static InputParameters validParams();

  MeltPoolHeight(const InputParameters & parameters);

  void initialize() override;
  void execute() override;
  void threadJoin(const UserObject & uo) override;
  virtual Real getValue() override;

protected:
  virtual void computeQpValue();
  /// Level set variable
  const VariableValue & _ls;
  const MaterialProperty<Real> & _f_l;
  Real _min_value;
  Real _max_value;
  Real _height;
  /// Current quadrature point
  unsigned int _qp;
};
