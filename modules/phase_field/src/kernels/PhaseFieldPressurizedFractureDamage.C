//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PhaseFieldPressurizedFractureDamage.h"

registerMooseObject("PhaseFieldApp", PhaseFieldPressurizedFractureDamage);

InputParameters
PhaseFieldPressurizedFractureDamage::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Kernel of phase field fracture damage variable due to pressure.");
  params.addRequiredCoupledVar("displacements",
                               "The string of displacements suitable for the problem statement");
  params.addParam<MaterialPropertyName>("mob_name", "L", "The mobility used with the kernel");
  return params;
}

PhaseFieldPressurizedFractureDamage::PhaseFieldPressurizedFractureDamage(
    const InputParameters & parameters)
  : DerivativeMaterialInterface<Kernel>(parameters),
    _ndisp(coupledComponents("displacements")),
    _grad_disp(3),
    _disp(3),
    _disp_var(3),
    _pressure(getDefaultMaterialProperty<Real>("fracture_pressure")),
    _L(getMaterialProperty<Real>("mob_name"))
{
  for (unsigned int i = 0; i < _ndisp; ++i)
  {
    _grad_disp[i] = &coupledGradient("displacements", i);
    _disp[i] = &coupledValue("displacements", i);
    _disp_var[i] = coupled("displacements", i);
  }
  for (unsigned i = _ndisp; i < 3; ++i)
  {
    _disp[i] = &_zero;
    _grad_disp[i] = &_grad_zero;
  }
}

Real
PhaseFieldPressurizedFractureDamage::computeQpResidual()
{
  // return -_pressure[_qp] *
  //        ((*_grad_disp[0])[_qp](0) + (*_grad_disp[1])[_qp](1) + (*_grad_disp[2])[_qp](2)) *
  //        _test[_i][_qp];

  RealVectorValue u((*_disp[0])[_qp], (*_disp[1])[_qp], (*_disp[2])[_qp]);
  return _pressure[_qp] * u * _grad_test[_i][_qp] * _L[_qp];
}

Real
PhaseFieldPressurizedFractureDamage::computeQpJacobian()
{
  return 0.0;
}

Real
PhaseFieldPressurizedFractureDamage::computeQpOffDiagJacobian(unsigned int jvar)
{
  for (unsigned int coupled_component = 0; coupled_component < _ndisp; ++coupled_component)
    if (jvar == _disp_var[coupled_component])
      // return -_pressure[_qp] * _grad_test[_j][_qp](coupled_component) * _test[_i][_qp];
      return _pressure[_qp] * _phi[_j][_qp] * _grad_test[_i][_qp](coupled_component) * _L[_qp];

  return 0.0;
}
