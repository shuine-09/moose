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

class CrystalPlasticityGBMarkerUO;

template <>
InputParameters validParams<CrystalPlasticityGBMarkerUO>();

class CrystalPlasticityGBMarkerUO : public ElementUserObject
{
public:
  CrystalPlasticityGBMarkerUO(const InputParameters & parameters);

  const std::map<dof_id_type, Real> & getGBElementsMap() const { return _gb_elem_map; };

  void initialize() override{};
  void execute() override;
  void threadJoin(const UserObject & /*uo*/) override{};
  void finalize() override;

protected:
  std::map<dof_id_type, Real> _gb_elem_map;
  const MaterialProperty<RealVectorValue> & _Euler_angles_mat_prop;
};
