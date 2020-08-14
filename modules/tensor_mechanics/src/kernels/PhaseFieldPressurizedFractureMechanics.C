//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PhaseFieldPressurizedFractureMechanics.h"

registerMooseObject("TensorMechanicsApp", PhaseFieldPressurizedFractureMechanics);

InputParameters
PhaseFieldPressurizedFractureMechanics::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription(
      "Kernel of phase field fracture displacement variable due to pressure.");
  params.addCoupledVar("c", "Phase field damage variable.");
  params.addRequiredParam<unsigned int>("component",
                                        "An integer corresponding to the direction "
                                        "the variable this kernel acts in. (0 for x, "
                                        "1 for y, 2 for z)");
  return params;
}

PhaseFieldPressurizedFractureMechanics::PhaseFieldPressurizedFractureMechanics(
    const InputParameters & parameters)
  : DerivativeMaterialInterface<Kernel>(parameters),
    _component(getParam<unsigned int>("component")),
    _grad_c(coupledGradient("c")),
    _c_var(coupled("c")),
    _pressure(getDefaultMaterialProperty<Real>("fracture_pressure"))
{
}

Real
PhaseFieldPressurizedFractureMechanics::computeQpResidual()
{
  return _pressure[_qp] * _grad_c[_qp](_component) * _test[_i][_qp];
}

Real
PhaseFieldPressurizedFractureMechanics::computeQpJacobian()
{
  return 0.0;
}

Real
PhaseFieldPressurizedFractureMechanics::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _c_var)
    return _pressure[_qp] * _grad_test[_j][_qp](_component) * _test[_i][_qp];

  return 0.0;
}
