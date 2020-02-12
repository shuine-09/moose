//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DEDMass.h"

registerADMooseObject("NavierStokesApp", DEDMass);

defineADValidParams(
    DEDMass,
    ADKernelValue,
    params.addClassDescription("This class computes the mass equation residual and Jacobian "
                               "contributions (the latter using automatic differentiation) for the "
                               "incompressible Navier-Stokes "
                               "equations.");
    params.addRequiredParam<MaterialPropertyName>("rho", "Denisty");
    params.addRequiredParam<MaterialPropertyName>("mass_change", "Mass change");
    params.addRequiredCoupledVar("grad_level_set", "Regularized gradient of Level set variable"););

template <ComputeStage compute_stage>
DEDMass<compute_stage>::DEDMass(const InputParameters & parameters)
  : ADKernelValue<compute_stage>(parameters),
    _mass_strong_residual(getADMaterialProperty<Real>("mass_strong_residual")),
    _rho(getADMaterialProperty<Real>("rho")),
    _grad_rho(getADMaterialProperty<RealVectorValue>("grad_rho")),
    _mass_change(getADMaterialProperty<Real>("mass_change")),
    _powder_feed(getADMaterialProperty<Real>("powder_feed")),
    _grad_ls(adCoupledVectorValue("grad_level_set"))
{
}

template <ComputeStage compute_stage>
ADReal
DEDMass<compute_stage>::precomputeQpResidual()
{
  ADRealVectorValue normal = _grad_ls[_qp] / (_grad_ls[_qp] + RealVectorValue(1e-10)).norm();
  // return _mass_strong_residual[_qp] + _powder_feed[_qp] * normal * _grad_rho[_qp] / _rho[_qp];
  return _mass_strong_residual[_qp] +
         _powder_feed[_qp] * (_grad_ls[_qp] + RealVectorValue(1e-10)).norm() / _rho[_qp];
  // return _mass_strong_residual[_qp] + 2 * _mass_change[_qp] / _rho[_qp];
}
