//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VacancyFluxZrPostprocessorForOpti.h"

registerMooseObject("XFEMApp", VacancyFluxZrPostprocessorForOpti);

InputParameters
VacancyFluxZrPostprocessorForOpti::validParams()
{
  InputParameters params = GeneralPostprocessor::validParams();
  params.addClassDescription(
      "Retrieves the oxygen vacancy flux Jv in the oxide.");
  params.addRequiredParam<UserObjectName>("velocity_uo",
      "The name of the userobject that obtains the velocity of the oxide interface and computes Jv.");
  return params;
}

VacancyFluxZrPostprocessorForOpti::VacancyFluxZrPostprocessorForOpti(const InputParameters & parameters)
  : GeneralPostprocessor(parameters)
{
}

void
VacancyFluxZrPostprocessorForOpti::initialize()
{
  const UserObject * uo =
      &(_fe_problem.getUserObjectBase(getParam<UserObjectName>("velocity_uo")));

  if (dynamic_cast<const XFEMC4VelocityZrOxAForOpti *>(uo) == nullptr)
    mooseError("UserObject casting to XFEMC4VelocityZrOxA in VacancyFluxZrPostprocessor");

  _velocity_uo = dynamic_cast<const XFEMC4VelocityZrOxAForOpti *>(uo);
}


Real
VacancyFluxZrPostprocessorForOpti::getValue()
{
  return _velocity_uo->getVacancyFlux();
}
