//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LevelSetTriMaterialReal.h"

registerMooseObject("XFEMApp", LevelSetTriMaterialReal);

InputParameters
LevelSetTriMaterialReal::validParams()
{
  InputParameters params = LevelSetTriMaterialBase::validParams();
  params.addClassDescription(
      "Compute a Real material property for tri-materials problem (consisting of three "
      "different materials) defined by 2 level set functions.");
  return params;
}

LevelSetTriMaterialReal::LevelSetTriMaterialReal(const InputParameters & parameters)
  : LevelSetTriMaterialBase(parameters),
    _bimaterial_material_prop(3),
    _material_prop(declareProperty<Real>(_base_name + _prop_name))
{
  _bimaterial_material_prop[0] = &getMaterialProperty<Real>(
      getParam<std::string>("levelset_neg_neg_base") + "_" + _prop_name);
  _bimaterial_material_prop[1] = &getMaterialProperty<Real>(
      getParam<std::string>("levelset_pos_neg_base") + "_" + _prop_name);
  _bimaterial_material_prop[2] = &getMaterialProperty<Real>(
      getParam<std::string>("levelset_pos_pos_base") + "_" + _prop_name);
}

void
LevelSetTriMaterialReal::assignQpPropertiesForLevelSetNegNeg()
{
  _material_prop[_qp] = (*_bimaterial_material_prop[0])[_qp];
}

void
LevelSetTriMaterialReal::assignQpPropertiesForLevelSetPosNeg()
{
  _material_prop[_qp] = (*_bimaterial_material_prop[1])[_qp];
}

void
LevelSetTriMaterialReal::assignQpPropertiesForLevelSetPosPos()
{
  _material_prop[_qp] = (*_bimaterial_material_prop[2])[_qp];
}
