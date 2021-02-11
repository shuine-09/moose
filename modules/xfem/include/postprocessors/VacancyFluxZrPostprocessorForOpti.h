#pragma once

#include "XFEMC4VelocityZrOxAForOpti.h"
#include "GeneralPostprocessor.h"

/**
*
* Retrieves the vacancy flux calculated by the C4 model for Zr corrosion
*
*/


class VacancyFluxZrPostprocessorForOpti : public GeneralPostprocessor
{
public:
  static InputParameters validParams();

  VacancyFluxZrPostprocessorForOpti(const InputParameters & parameters);

  virtual void initialize() override;

  virtual void execute() override {}

  virtual Real getValue() override;

protected:
  /// Pointer to XFEMC4VelocityZrOxA object
  const XFEMC4VelocityZrOxAForOpti * _velocity_uo;

};
