#pragma once

#include "ElementPostprocessor.h"

/**
 * Computes the average interface velocity.
 */
class AverageInterfaceVelocity : public ElementPostprocessor
{
public:
  static InputParameters validParams();

  AverageInterfaceVelocity(const InputParameters & parameters);
  void initialize() override;
  void execute() override;
  void finalize() override;
  void threadJoin(const UserObject & user_object) override;
  virtual PostprocessorValue getValue() override;

private:
  /// The max velocity on an element, this is done simply to avoid creating temporary calls to execute.
  Real _max_velocity;

  /// Level set delta function
  const ADMaterialProperty<Real> & _delta_function;

  /// Velocity vector variable
  const ADVectorVariableValue & _velocity;
};
