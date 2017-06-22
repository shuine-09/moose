/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "PFFracBulkRate.h"
#include "MathUtils.h"

template <>
InputParameters
validParams<PFFracBulkRate>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription(
      "Kernel to compute bulk energy contribution to damage order parameter residual equation");
  params.addRequiredParam<Real>("l", "Interface width");
  params.addRequiredParam<MaterialPropertyName>("gc_prop_var",
                                                "Material property name with gc value");
  params.addRequiredParam<MaterialPropertyName>(
      "G0_var", "Material property name with undamaged strain energy driving damage (G0_pos)");
  params.addParam<MaterialPropertyName>(
      "dG0_dstrain_var", "Material property name with derivative of G0_pos with strain");
  params.addCoupledVar("displacements",
                       "The string of displacements suitable for the problem statement");
  params.addParam<std::string>("base_name", "Material property base name");

  return params;
}

PFFracBulkRate::PFFracBulkRate(const InputParameters & parameters)
  : Kernel(parameters),
    _gc_prop(getMaterialProperty<Real>("gc_prop_var")),
    _G0_pos(getMaterialProperty<Real>("G0_var")),
    _dG0_pos_dstrain(isParamValid("dG0_dstrain_var")
                         ? &getMaterialProperty<RankTwoTensor>("dG0_dstrain_var")
                         : NULL),
    _ndisp(coupledComponents("displacements")),
    _disp_var(_ndisp),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _l(getParam<Real>("l"))
{
  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_var[i] = coupled("displacements", i);
}

Real
PFFracBulkRate::computeQpResidual()
{
  return -_gc_prop[_qp] * _l * _grad_u[_qp] * _grad_test[_i][_qp] +
         2.0 * (1.0 - _u[_qp]) * _test[_i][_qp] * _G0_pos[_qp] -
         _gc_prop[_qp] / _l * _u[_qp] * _test[_i][_qp];
}

Real
PFFracBulkRate::computeQpJacobian()
{
  return -_gc_prop[_qp] * _l * _grad_phi[_j][_qp] * _grad_test[_i][_qp] -
         2.0 * _test[_j][_qp] * _test[_i][_qp] * _G0_pos[_qp] -
         _gc_prop[_qp] / _l * _phi[_j][_qp] * _test[_i][_qp];
}
