#pragma once

#include "ADKernelGrad.h"

/**
 * Implements the re-initialization equation that uses regularized gradient.
 */
class LevelSetGradientRegularizationReinitializationSUPG : public ADKernelGrad
{
public:
  static InputParameters validParams();

  LevelSetGradientRegularizationReinitializationSUPG(const InputParameters & parameters);

protected:
  virtual ADRealVectorValue precomputeQpResidual() override;

  /// Regularized gradient of the level set variable at time, \tau = 0.
  const ADVectorVariableValue & _grad_c;

  /// Interface thickness
  const Real & _epsilon;

  /// factor
  const Real & _gamma;

  /// Velocity vector variable
  const ADVectorVariableValue & _velocity;
};
