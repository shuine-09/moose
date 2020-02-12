//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Diffusion.h"

registerMooseObject("MooseApp", Diffusion);

defineLegacyParams(Diffusion);

InputParameters
Diffusion::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("The Laplacian operator ($-\\nabla \\cdot \\nabla u$), with the weak "
                             "form of $(\\nabla \\phi_i, \\nabla u_h)$.");
  params.addRequiredCoupledVar(
      "displacements",
      "The displacements appropriate for the simulation geometry and coordinate system");
  return params;
}

Diffusion::Diffusion(const InputParameters & parameters)
  : Kernel(parameters), _ndisp(coupledComponents("displacements")), _disp(3)
{
  for (unsigned int i = 0; i < _ndisp; ++i)
  {
    _disp[i] = &coupledValue("displacements", i);
  }
  // set unused dimensions to zero
  for (unsigned i = _ndisp; i < 3; ++i)
  {
    _disp[i] = &_zero;
  }
}

Real
Diffusion::computeQpResidual()
{
  RealVectorValue disp((*_disp[0])[_qp], (*_disp[1])[_qp], (*_disp[2])[_qp]);
  return _grad_u[_qp] * disp * (1.0e-3);
}

Real
Diffusion::computeQpJacobian()
{
  return 0.0;
}
