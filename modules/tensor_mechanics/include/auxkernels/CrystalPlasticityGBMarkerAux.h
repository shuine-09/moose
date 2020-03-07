//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"
#include "CrystalPlasticityGBMarkerUO.h"

class CrystalPlasticityGBMarkerAux;

template <>
InputParameters validParams<CrystalPlasticityGBMarkerAux>();

class CrystalPlasticityGBMarkerAux : public AuxKernel
{
public:
  /**
   * CrystalPlasticityGBMarkerAux
   */
  CrystalPlasticityGBMarkerAux(const InputParameters & parameters);

  virtual ~CrystalPlasticityGBMarkerAux() {}

protected:
  virtual Real computeValue();

  const CrystalPlasticityGBMarkerUO * _marker_uo;
  const std::map<dof_id_type, Real> * _marker_map;
  bool _second_layer;
};
