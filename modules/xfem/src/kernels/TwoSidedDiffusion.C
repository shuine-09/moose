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

#include "TwoSidedDiffusion.h"
#include "FEProblem.h"
#include "DisplacedProblem.h"
#include "Assembly.h"
#include "XFEM.h"
#include "MooseMesh.h"

template<>
InputParameters validParams<TwoSidedDiffusion>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("The Laplacian operator ($-\\nabla \\cdot \\nabla u$), with the weak form of $(\\nabla \\phi_i, \\nabla u_h)$.");
  params.addParam<Real>("diffusivity_ls_plus", 1.0, "Diffusivity Coefficient for domain of positive level set");
  params.addParam<Real>("diffusivity_ls_minus", 1.0, "Diffusivity Coefficient for domain of minus level set");
  params.addRequiredCoupledVar("level_set", "Level set variable");
  return params;
}

TwoSidedDiffusion::TwoSidedDiffusion(const InputParameters & parameters) :
    Kernel(parameters),
    _diffusivity_ls_plus(getParam<Real>("diffusivity_ls_plus")),
    _diffusivity_ls_minus(getParam<Real>("diffusivity_ls_minus")),
    _ls(coupledValue("level_set"))
{
  FEProblem * fe_problem = dynamic_cast<FEProblem *>(&_subproblem);
  if (fe_problem == NULL)
    mooseError("Problem casting _subproblem to FEProblem in TwoSidedDiffusion");
  _xfem = MooseSharedNamespace::dynamic_pointer_cast<XFEM>(fe_problem->getXFEM());
  if (_xfem == NULL)
    mooseError("Problem casting to XFEM in TwoSidedDiffusion");
}

bool 
TwoSidedDiffusion::isPlusDomain()
{
  bool isPlus = true;
  const Elem * undisplaced_elem  = NULL;
  FEProblem * _fe_problem = dynamic_cast<FEProblem *>(&_subproblem);
  if(_fe_problem->getDisplacedProblem() != NULL)
    undisplaced_elem = _fe_problem->getDisplacedProblem()->refMesh().elemPtr(_current_elem->id());
  else
    undisplaced_elem = _current_elem;

  for ( unsigned int qp = 0; qp < _q_point.size(); ++qp)
  {
    Real flag = _xfem->flagQpointInside(undisplaced_elem, _q_point[qp]); //inside (flag = 1) or ouside (flag = 0) real domain
    if (flag > 0.5)
    {
      if (_ls[qp] > 0.0)
        isPlus = false;
      else
        isPlus = true;

      break;
    }
  }
  return isPlus;
}

Real
TwoSidedDiffusion::computeQpResidual()
{
  Real diffusivity = 0.0;

  if (isPlusDomain())
    diffusivity = _diffusivity_ls_plus;
  else
    diffusivity = _diffusivity_ls_minus;

  return _grad_u[_qp] * _grad_test[_i][_qp] * diffusivity;
}

Real
TwoSidedDiffusion::computeQpJacobian()
{
  Real diffusivity = 0.0;

  if (isPlusDomain())
    diffusivity = _diffusivity_ls_plus;
  else
    diffusivity = _diffusivity_ls_minus;

  return _grad_phi[_j][_qp] * _grad_test[_i][_qp] * diffusivity;
}
