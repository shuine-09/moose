#include "LevelSetGradientRegularizationReinitializationSUPG.h"

registerMooseObject("NavierStokesApp", LevelSetGradientRegularizationReinitializationSUPG);

InputParameters
LevelSetGradientRegularizationReinitializationSUPG::validParams()
{
  InputParameters params = ADKernelGrad::validParams();
  params.addClassDescription("The re-initialization equation that uses regularized gradient.");
  params.addRequiredCoupledVar("level_set_gradient",
                               "Regularized gradient of the level set variable");
  params.addRequiredParam<Real>(
      "epsilon", "The epsilon coefficient to be used in the reinitialization calculation.");
  params.addParam<PostprocessorName>("gamma", 1, "Reinitilzation parameter.");
  params.addRequiredCoupledVar("velocity", "Velocity vector variable.");
  return params;
}

LevelSetGradientRegularizationReinitializationSUPG::
    LevelSetGradientRegularizationReinitializationSUPG(const InputParameters & parameters)
  : ADKernelGrad(parameters),
    _grad_c(adCoupledVectorValue("level_set_gradient")),
    _epsilon(getParam<Real>("epsilon")),
    _gamma(getPostprocessorValue("gamma")),
    _velocity(adCoupledVectorValue("velocity"))
{
}

ADRealVectorValue
LevelSetGradientRegularizationReinitializationSUPG::precomputeQpResidual()
{
  ADReal s = (_grad_c[_qp] + ADRealVectorValue(libMesh::TOLERANCE)).norm();
  ADRealVectorValue n_hat = _grad_c[_qp] / s;
  ADRealVectorValue f = _u[_qp] * (1 - _u[_qp]) * n_hat;
  ADReal tau =
      _current_elem->hmin() /
      (2 * (_velocity[_qp] + RealVectorValue(libMesh::TOLERANCE * libMesh::TOLERANCE)).norm());
  return (tau * _velocity[_qp]) * _gamma * (-f + _epsilon * _grad_u[_qp]);
}
