/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "XFEMWeibullAux.h"

#include "XFEM.h"

template<>
InputParameters validParams<XFEMWeibullAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<Real>("weibull_modulus","The Weibull modulus quantifying observed variability in the property.");
  params.addRequiredParam<Real>("specimen_material_property", "The median value of material property observed in the laboratory.");
  params.addRequiredParam<Real>("specimen_volume", "Specimen volume used in the laboratory.");
  return params;
}

XFEMWeibullAux::XFEMWeibullAux(const InputParameters & parameters)
  :AuxKernel(parameters),
  _weibull_modulus(getParam<Real>("weibull_modulus")),
  _specimen_volume(getParam<Real>("specimen_volume")),
  _specimen_material_property(getParam<Real>("specimen_material_property"))
{
  FEProblem * fe_problem = dynamic_cast<FEProblem *>(&_subproblem);
  if (fe_problem == NULL)
    mooseError("Problem casting _subproblem to FEProblem in XFEMWeibullAux");
  _xfem = fe_problem->get_xfem();
  if (isNodal())
    mooseError("XFEMWeibullAux can only be run on an element variable");

  // Setup the random number generation
  setRandomResetFrequency(EXEC_TIMESTEP_BEGIN);
}

Real
XFEMWeibullAux::computeValue()
{
  Real rn = getRandomReal();
  Real eta = -1.0;
  if ( _weibull_modulus != 0)
    eta = _specimen_material_property * std::pow(_specimen_volume*std::log(rn)/(_current_elem_volume*std::log(0.5)),1.0/(Real)_weibull_modulus);
  return eta;
}
