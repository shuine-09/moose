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
#include "XFEMWeibullRandom.h"
#include "XFEM.h"

template<>
InputParameters validParams<XFEMWeibullRandom>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<Real>("weibull_modulus","The Weibull modulus quantifying observed variability in the property.");
  params.addRequiredParam<Real>("specimen_material_property", "The median value of material property observed in the laboratory.");
  params.addRequiredParam<Real>("specimen_volume", "Specimen volume used in the laboratory.");
  return params;
}

XFEMWeibullRandom::XFEMWeibullRandom(const InputParameters & parameters)
  :Material(parameters),
  _weibull_eta(declareProperty<Real>("weibull_eta")),
  _weibull_eta_old(declarePropertyOld<Real>("weibull_eta")),
  _weibull_modulus(getParam<Real>("weibull_modulus")),
  _specimen_volume(getParam<Real>("specimen_volume")),
  _specimen_material_property(getParam<Real>("specimen_material_property"))
{
  FEProblem * fe_problem = dynamic_cast<FEProblem *>(&_subproblem);

  if (fe_problem == NULL)
    mooseError("Problem casting _subproblem to FEProblem in XFEMWeibullRandom");
  _xfem = MooseSharedNamespace::dynamic_pointer_cast<XFEM>(fe_problem->getXFEM());
  if (_xfem == NULL)
    mooseError("Problem casting to XFEM in XFEMWeibullRandom");

  // Setup the random number generation
  setRandomResetFrequency(EXEC_INITIAL);
}

void
XFEMWeibullRandom::initQpStatefulProperties()
{
  if (_qp == 0)
  {
    if ( std::abs(_weibull_modulus) > 1.0e-5)
    {
      Real rn = getRandomReal();
      _eta = _specimen_material_property * std::pow(_specimen_volume*std::log(rn)/(_current_elem->volume()*std::log(0.5)),1.0/(Real)_weibull_modulus);
    }
    _weibull_eta[_qp] = _eta;
    _weibull_eta_old[_qp] = _eta;
  }
  else
  {
    _weibull_eta[_qp] = _eta;
    _weibull_eta_old[_qp] = _eta;
  }
}

void
XFEMWeibullRandom::computeQpProperties()
{
}
