//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
/**
#include "C4VacancyFluxPostprocessor.h"

registerMooseObject("XFEMApp", C4VacancyFluxPostprocessor);

InputParameters
C4VacancyFluxPostprocessor::validParams()
{
  InputParameters params = GeneralPostprocessor::validParams();
  params.addClassDescription(
      "Retrieves the value of the vacancy flux of the C4 model ");
  params.addRequiredParam<UserObjectName>("velocity_uo",
      "The name of the userobject that obtains the velocity of the oxide/metal interface.");
  return params;
}

C4VacancyFluxPostprocessor::C4VacancyFluxPostprocessor(const InputParameters & parameters)
  : GeneralPostprocessor(parameters)
{
}

void
C4VacancyFluxPostprocessor::initialize()
{
  const UserObject * uo =
      &(_fe_problem.getUserObjectBase(getParam<UserObjectName>("velocity_uo")));

  if (dynamic_cast<const XFEMC4VelocityOxideWeakMicro *>(uo) == nullptr)
    mooseError("UserObject casting to XFEMC4VelocityOxideWeakMicro in C4VacancyFluxPostprocessor");

  _velocity_uo = dynamic_cast<const XFEMC4VelocityOxideWeakMicro *>(uo);
}


Real
C4VacancyFluxPostprocessor::getValue()
{
  return _velocity_uo->getVacancyFlux();
}
*/
