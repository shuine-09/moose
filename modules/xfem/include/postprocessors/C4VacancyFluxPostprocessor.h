/**
#pragma once

#include "XFEMC4VelocityOxideWeakMicro.h"
#include "GeneralPostprocessor.h"

class C4VacancyFluxPostprocessor : public GeneralPostprocessor
{
public:
  static InputParameters validParams();

  C4VacancyFluxPostprocessor(const InputParameters & parameters);

  virtual void initialize() override;

  virtual void execute() override {}

  virtual Real getValue() override;

protected:
  /// Pointer to PointValueAtXFEMInterface object
  const XFEMC4VelocityOxideWeakMicro * _velocity_uo;

};
*/
