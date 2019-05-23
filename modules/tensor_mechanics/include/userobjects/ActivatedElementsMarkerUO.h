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

class ActivatedElementsMarkerUO;

template <>
InputParameters validParams<ActivatedElementsMarkerUO>();

class ActivatedElementsMarkerUO : public ElementUserObject
{
public:
  ActivatedElementsMarkerUO(const InputParameters & parameters);

  const std::map<dof_id_type, Real> & getActivatedElementsMap() const
  {
    return _activated_elem_map;
  };

  void initialize() override{};
  void execute() override;
  void threadJoin(const UserObject & /*uo*/) override{};
  void finalize() override;

protected:
  std::map<dof_id_type, Real> _activated_elem_map;
  const VariableValue & _temp_aux;
  Real _melt_temperature;
};
