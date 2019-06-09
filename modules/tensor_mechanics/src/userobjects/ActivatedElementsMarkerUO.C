//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ActivatedElementsMarkerUO.h"
#include "libmesh/quadrature.h"
#include "libmesh/parallel_algebra.h"
#include "libmesh/parallel.h"

registerMooseObject("MooseApp", ActivatedElementsMarkerUO);

template <>
InputParameters
validParams<ActivatedElementsMarkerUO>()
{
  InputParameters params = validParams<ElementUserObject>();
  params.addClassDescription("Determine activated elements.");
  params.addRequiredCoupledVar("temp_aux",
                               "Temperature aux variable used to determine activated elements.");
  params.addRequiredParam<Real>("melt_temperature", "Melt temperature.");

  return params;
}

ActivatedElementsMarkerUO::ActivatedElementsMarkerUO(const InputParameters & parameters)
  : ElementUserObject(parameters),
    _temp_aux(coupledValue("temp_aux")),
    _melt_temperature(getParam<Real>("melt_temperature"))
{
}

void
ActivatedElementsMarkerUO::execute()
{
  Real marker_old = 0;
  auto marker = _activated_elem_map.find(_current_elem->id());
  if (marker != _activated_elem_map.end())
    marker_old = marker->second;

  Real temp = 0;
  for (unsigned int qp = 0; qp < _qrule->n_points(); ++qp)
    temp += _temp_aux[qp];

  temp /= _qrule->n_points();

  if (temp > _melt_temperature || marker_old >= 1.0 || _current_elem->subdomain_id() != 1)
    _activated_elem_map[_current_elem->id()] = 1.0;
  else
    _activated_elem_map[_current_elem->id()] = 0.0;
}

void
ActivatedElementsMarkerUO::finalize()
{
  _communicator.set_union(_activated_elem_map);
}
