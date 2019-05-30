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
#include "ActivatedElementsMarkerUO.h"

class ActivatedElementsMarker;

template <>
InputParameters validParams<ActivatedElementsMarker>();

class ActivatedElementsMarker : public AuxKernel
{
public:
  /**
   * ActivatedElementsMarker
   */
  ActivatedElementsMarker(const InputParameters & parameters);

  virtual ~ActivatedElementsMarker() {}

protected:
  virtual Real computeValue();

  const VariableValue & _temp_aux;
  Real _melt_temperature;

  // const ActivatedElementsMarkerUO * _marker_uo;
  // const std::map<dof_id_type, Real> * _marker_map;
};
