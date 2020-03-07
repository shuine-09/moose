//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CrystalPlasticityGBMarkerAux.h"

registerMooseObject("TensorMechanicsApp", CrystalPlasticityGBMarkerAux);

template <>
InputParameters
validParams<CrystalPlasticityGBMarkerAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addParam<UserObjectName>("marker_uo", "Marker UserObject");
  params.addParam<bool>("second_layer", false, "Add second layer.");
  return params;
}

CrystalPlasticityGBMarkerAux::CrystalPlasticityGBMarkerAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _marker_uo(isParamValid("marker_uo") ? &getUserObjectByName<CrystalPlasticityGBMarkerUO>(
                                               getParam<UserObjectName>("marker_uo"))
                                         : nullptr),
    _second_layer(getParam<bool>("second_layer"))

{
  if (_marker_uo)
    _marker_map = &(_marker_uo->getGBElementsMap());
  else
    _marker_map = nullptr;
}

Real
CrystalPlasticityGBMarkerAux::computeValue()
{
  dof_id_type elem_id = _current_elem->id();

  Real gb0 = _marker_map->find(elem_id)->second;
  bool isGB = false;

  if (!_second_layer)
  {
    for (unsigned int i = 0; i < _current_elem->n_neighbors(); ++i)
    {
      if (_current_elem->neighbor(i))
      {
        Real gbi = _marker_map->find(_current_elem->neighbor(i)->id())->second;
        if (gbi != gb0)
        {
          isGB = true;
          break;
        }
      }
    }
  }
  else
  {
    for (unsigned int i = 0; i < _current_elem->n_neighbors(); ++i)
    {
      if (_current_elem->neighbor(i))
      {
        Real gbi = _marker_map->find(_current_elem->neighbor(i)->id())->second;
        if (gbi != gb0)
        {
          isGB = true;
          break;
        }

        const Elem * neigh_elem = _current_elem->neighbor(i);
        for (unsigned int j = 0; j < neigh_elem->n_neighbors(); ++j)
        {
          if (neigh_elem->neighbor(j))
          {
            Real gbj = _marker_map->find(neigh_elem->neighbor(j)->id())->second;
            if (gbj != gb0)
            {
              isGB = true;
              break;
            }
          }
        }
      }
    }
  }

  if (isGB)
    return 1;
  else
    return 0;
}
