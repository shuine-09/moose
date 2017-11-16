/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "LevelSetMultiStressMaterial.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

template <>
InputParameters
validParams<LevelSetMultiStressMaterial>()
{
  InputParameters params = validParams<Material>();
  params.addClassDescription("Compute a global stress form multiple phase stresses");
  params.addRequiredCoupledVar("level_set_var", "Level set variable");
  params.addRequiredParam<std::string>("levelset_plus_base",
                                       "Base name for the strain in level set plus region.");
  params.addRequiredParam<std::string>("levelset_minus_base",
                                       "Base name for the strain in level set minus region.");
  params.addParam<std::string>("base_name", "Base name for the computed global stress (optional)");
  return params;
}

LevelSetMultiStressMaterial::LevelSetMultiStressMaterial(const InputParameters & parameters)
  : Material(parameters),
    _ls(coupledValue("level_set_var")),
    _phase_base(2),
    _phase_stress(2),
    _dphase_stress_dstrain(2),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _stress(declareProperty<RankTwoTensor>(_base_name + "stress")),
    _dstress_dstrain(declareProperty<RankFourTensor>(_base_name + "Jacobian_mult"))
{
  _phase_stress[0] =
      &getMaterialProperty<RankTwoTensor>(getParam<std::string>("levelset_plus_base") + "_stress");
  _dphase_stress_dstrain[0] = &getMaterialProperty<RankFourTensor>(
      getParam<std::string>("levelset_plus_base") + "_Jacobian_mult");

  _phase_stress[1] =
      &getMaterialProperty<RankTwoTensor>(getParam<std::string>("levelset_minus_base") + "_stress");
  _dphase_stress_dstrain[1] = &getMaterialProperty<RankFourTensor>(
      getParam<std::string>("levelset_minus_base") + "_Jacobian_mult");
}

void
LevelSetMultiStressMaterial::computeQpProperties()
{
  _stress[_qp].zero();
  _dstress_dstrain[_qp].zero();

  if (_ls[_qp] > 0)
  {
    _stress[_qp] = (*_phase_stress[0])[_qp];
    _dstress_dstrain[_qp] = (*_dphase_stress_dstrain[0])[_qp];
  }
  else
  {
    _stress[_qp] = (*_phase_stress[1])[_qp];
    _dstress_dstrain[_qp] = (*_dphase_stress_dstrain[1])[_qp];
  }
}
