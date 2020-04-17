//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MeltPoolLevelSetINSMaterial.h"
#include "Function.h"

registerADMooseObject("NavierStokesApp", MeltPoolLevelSetINSMaterial);

defineADValidParams(
    MeltPoolLevelSetINSMaterial,
    INSADTauMaterial,
    params.addRequiredCoupledVar("level_set", "Level set variable");
    params.addRequiredCoupledVar("grad_level_set", "Regularized gradient of Level set variable");
    params.addRequiredCoupledVar("temperature", "Temperature variable");
    params.addRequiredCoupledVar("curvature", "Level set variable");
    params.addRequiredParam<Real>("surface_tension", "surface tension coefficient.");
    params.addRequiredParam<Real>("thermal_expansion", "thermal expansion.");
    params.addParam<Real>("reference_temperature", 300, "reference_temperature.");
    params.addRequiredParam<Real>("thermal_capillary", "thermalcapillary coefficient.");
    params.addParam<UserObjectName>("location", "User object name of laser location");
    params.addRequiredParam<Real>("rho_l", "Liquid density."););

template <ComputeStage compute_stage>
MeltPoolLevelSetINSMaterial<compute_stage>::MeltPoolLevelSetINSMaterial(
    const InputParameters & parameters)
  : INSADTauMaterial<compute_stage>(parameters),
    _ls(adCoupledValue("level_set")),
    _grad_ls(adCoupledVectorValue("grad_level_set")),
    _temp(adCoupledValue("temperature")),
    _grad_temp(adCoupledGradient("temperature")),
    _curvature(adCoupledValue("curvature")),
    _thermal_expansion(getParam<Real>("thermal_expansion")),
    _reference_temperature(getParam<Real>("reference_temperature")),
    _permeability(getADMaterialProperty<Real>("permeability")),
    _sigma(getParam<Real>("surface_tension")),
    _sigmaT(getParam<Real>("thermal_capillary")),
    _rho_l(getParam<Real>("rho_l")),
    _melt_pool_momentum_source(declareADProperty<RealVectorValue>("melt_pool_momentum_source")),
    _location(getUserObjectByName<MeltPoolLevelSetLocation>("location"))
{
}

template <ComputeStage compute_stage>
void
MeltPoolLevelSetINSMaterial<compute_stage>::computeQpProperties()
{
  INSADTauMaterial<compute_stage>::computeQpProperties();

  ADRealVectorValue normal = _grad_ls[_qp] / (_grad_ls[_qp] + RealVectorValue(1e-10)).norm();
  RankTwoTensor iden(RankTwoTensor::initIdentity);
  ADRankTwoTensor proj;
  proj.vectorOuterProduct(normal, normal);
  proj = iden - proj;

  ADRealVectorValue darcy_term = -_permeability[_qp] * std::abs(1 - _ls[_qp]) * _velocity[_qp];
  ADRealVectorValue surface_tension_term = _sigma * _curvature[_qp] * _grad_ls[_qp];
  ADRealVectorValue thermalcapillary_term =
      -proj * _grad_temp[_qp] * _sigmaT * (_grad_ls[_qp] + RealVectorValue(1e-10)).norm();
  ADRealVectorValue bouyancy_term = -_rho_l * _gravity * _thermal_expansion *
                                    (_temp[_qp] - _reference_temperature) * (1 - _ls[_qp]);

  _melt_pool_momentum_source[_qp] =
      -thermalcapillary_term + surface_tension_term + bouyancy_term + darcy_term;

  _momentum_strong_residual[_qp] -= _melt_pool_momentum_source[_qp];
}
