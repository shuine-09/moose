/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "LevelSetBiMaterialProperty.h"

template <>
InputParameters
validParams<LevelSetBiMaterialProperty>()
{
  InputParameters params = validParams<LevelSetBiMaterialBase>();
  params.addClassDescription(
      "Compute the diffusion coefficient for two materials defined by a level set function.");
  params.addRequiredParam<std::string>("prop_name",
                                       "Name for the computed material property (optional)");
  return params;
}

LevelSetBiMaterialProperty::LevelSetBiMaterialProperty(const InputParameters & parameters)
  : LevelSetBiMaterialBase(parameters),
    _bimaterial_material_prop(2),
    _prop_name(getParam<std::string>("prop_name")),
    _material_prop(declareProperty<Real>(_base_name + _prop_name))
{
  _bimaterial_material_prop[0] = &getMaterialProperty<Real>(
      getParam<std::string>("levelset_positive_base") + "_" + _prop_name);
  _bimaterial_material_prop[1] = &getMaterialProperty<Real>(
      getParam<std::string>("levelset_negative_base") + "_" + _prop_name);
}

void
LevelSetBiMaterialProperty::assignQpPropertiesForLevelSetPositive()
{
  _material_prop[_qp] = (*_bimaterial_material_prop[0])[_qp];
}

void
LevelSetBiMaterialProperty::assignQpPropertiesForLevelSetNegative()
{
  _material_prop[_qp] = (*_bimaterial_material_prop[1])[_qp];
}
