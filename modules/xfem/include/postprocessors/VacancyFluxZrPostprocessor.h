#pragma once

#include "XFEMC4VelocityZrOxA.h"
#include "GeneralPostprocessor.h"

/**
*
* Retrieves the vacancy flux calculated by the C4 model for Zr corrosion
*
*/


class VacancyFluxZrPostprocessor : public GeneralPostprocessor
{
public:
  static InputParameters validParams();

  VacancyFluxZrPostprocessor(const InputParameters & parameters);

  virtual void initialize() override;

  virtual void execute() override {}

  virtual Real getValue() override;

protected:
  /// Pointer to XFEMC4VelocityZrOxA object
  const XFEMC4VelocityZrOxA * _velocity_uo;

};
