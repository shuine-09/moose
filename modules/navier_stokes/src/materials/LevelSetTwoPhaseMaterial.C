/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "LevelSetTwoPhaseMaterial.h"

registerADMooseObject("NavierStokesApp", LevelSetTwoPhaseMaterial);

defineADValidParams(
    LevelSetTwoPhaseMaterial, ADMaterial, params.addRequiredCoupledVar("phi", "Level set variable");
    params.addRequiredParam<MaterialPropertyName>("prop_name", "Name of material property.");
    params.addRequiredParam<Real>("prop_value_1", "value of phase 1.");
    params.addRequiredParam<Real>("prop_value_2", "value of phase 2."););

template <ComputeStage compute_stage>
LevelSetTwoPhaseMaterial<compute_stage>::LevelSetTwoPhaseMaterial(
    const InputParameters & parameters)
  : ADMaterial<compute_stage>(parameters),
    _phi(adCoupledValue("phi")),
    _prop_value_1(getParam<Real>("prop_value_1")),
    _prop_value_2(getParam<Real>("prop_value_2")),
    _prop(declareADProperty<Real>(getParam<MaterialPropertyName>("prop_name")))
{
}

template <ComputeStage compute_stage>
void
LevelSetTwoPhaseMaterial<compute_stage>::computeQpProperties()
{
  _prop[_qp] = _prop_value_1 + (_prop_value_2 - _prop_value_1) * _phi[_qp];
}
