//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CrystalPlasticityGBMarkerUO.h"
#include "libmesh/quadrature.h"
#include "libmesh/parallel_algebra.h"
#include "libmesh/parallel.h"

registerMooseObject("MooseApp", CrystalPlasticityGBMarkerUO);

template <>
InputParameters
validParams<CrystalPlasticityGBMarkerUO>()
{
  InputParameters params = validParams<ElementUserObject>();
  return params;
}

CrystalPlasticityGBMarkerUO::CrystalPlasticityGBMarkerUO(const InputParameters & parameters)
  : ElementUserObject(parameters),
    _Euler_angles_mat_prop(getMaterialProperty<RealVectorValue>("Euler_angles"))
{
}

void
CrystalPlasticityGBMarkerUO::execute()
{
  _gb_elem_map[_current_elem->id()] = _Euler_angles_mat_prop[0](0);
}

void
CrystalPlasticityGBMarkerUO::finalize()
{
  _communicator.set_union(_gb_elem_map);
}
