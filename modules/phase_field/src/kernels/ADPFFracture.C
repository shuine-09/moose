//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "ADPFFracture.h"

registerADMooseObject("PhaseFieldApp", ADPFFracture);

defineADValidParams(
    ADPFFracture,
    ADKernel,
    params.addClassDescription(
        "Kernel to compute bulk energy contribution to damage order parameter residual equation");
    params.addParam<MaterialPropertyName>("l_name", "l", "Interface width");
    params.addParam<MaterialPropertyName>("gc", "gc_prop", "Critical fracture energy density"););

template <ComputeStage compute_stage>
ADPFFracture<compute_stage>::ADPFFracture(const InputParameters & parameters)
  : ADKernel<compute_stage>(parameters),
    _gc_prop(adGetMaterialProperty<Real>("gc")),
    _l(adGetMaterialProperty<Real>("l_name")),
    _hist(adGetADMaterialProperty<Real>("hist")),
    _hist_old(adGetMaterialPropertyOld<Real>("hist"))
{
}

template <ComputeStage compute_stage>
ADResidual
ADPFFracture<compute_stage>::computeQpResidual()
{

  const ADReal gc0 = _gc_prop[_qp];

  ADReal gc = gc0;

  if (_grad_u[_qp](0) * _grad_u[_qp](0) + _grad_u[_qp](1) * _grad_u[_qp](1) > 1.0e-4)
  {
    ADReal cosA = (_grad_u[_qp](0) * _grad_u[_qp](0) - _grad_u[_qp](1) * _grad_u[_qp](1)) /
                  (_grad_u[_qp](0) * _grad_u[_qp](0) + _grad_u[_qp](1) * _grad_u[_qp](1) + 1.e-4);
    ADReal sinA = (2.0 * _grad_u[_qp](0) * _grad_u[_qp](1)) /
                  (_grad_u[_qp](0) * _grad_u[_qp](0) + _grad_u[_qp](1) * _grad_u[_qp](1) + 1.e-4);

    ADReal cc = cosA * 0.5 + sinA * 0.8660254;

    gc = gc0 * (1.0 + 0.4 * cc) * (1.0 + 0.4 * cc);
  }

  // std::cout << "grad_u = " << _grad_u[_qp](0) << ", " << _grad_u[_qp](1) << std::endl;
  // std::cout << "theta = " << theta << std::endl;
  // std::cout << "gc = " << gc << std::endl;

  return (-gc * _l[_qp] * _grad_u[_qp] * _grad_test[_i][_qp] +
          2.0 * (1.0 - _u[_qp]) * _test[_i][_qp] * _hist_old[_qp] -
          gc / _l[_qp] * _u[_qp] * _test[_i][_qp]) /
         gc;
}

template class ADPFFracture<RESIDUAL>;
template class ADPFFracture<JACOBIAN>;
