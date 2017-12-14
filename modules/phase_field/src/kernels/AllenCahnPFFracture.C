/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "AllenCahnPFFracture.h"
#include "MathUtils.h"
#include "RankTwoTensor.h"

template <>
InputParameters
validParams<AllenCahnPFFracture>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription(
      "Kernel to compute bulk energy contribution to damage order parameter residual equation");
  params.addParam<MaterialPropertyName>("l_name", "l", "Interface width");
  params.addParam<MaterialPropertyName>("visco_name", "visco", "Viscosity parameter");
  params.addParam<MaterialPropertyName>("gc", "gc_prop", "Critical fracture energy density");
  // params.addRequiredCoupledVar("beta", "Variable storing the laplacian of c");
  params.addCoupledVar("displacements",
                       "The string of displacements suitable for the problem statement");
  params.addParam<MaterialPropertyName>(
      "F_name", "E_el", "Name of material property storing the elastic energy");

  return params;
}

AllenCahnPFFracture::AllenCahnPFFracture(const InputParameters & parameters)
  : DerivativeMaterialInterface<JvarMapKernelInterface<Kernel>>(parameters),
    _gc_prop(getMaterialProperty<Real>("gc")),
    // _beta(coupledValue("beta")),
    // _beta_var(coupled("beta")),
    _ndisp(coupledComponents("displacements")),
    _disp_var(_ndisp),
    _l(getMaterialProperty<Real>("l_name")),
    _visco(getMaterialProperty<Real>("visco_name")),
    _dFdc(getMaterialPropertyDerivative<Real>("F_name", _var.name())),
    _d2Fdc2(getMaterialPropertyDerivative<Real>("F_name", _var.name(), _var.name())),
    _d2Fdcdstrain(getMaterialProperty<RankTwoTensor>("d2Fdcdstrain")),
    _H0_pos(getMaterialProperty<Real>("H0_pos")),
    _H0_pos_old(getMaterialPropertyOld<Real>("H0_pos"))
{
  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_var[i] = coupled("displacements", i);
}

Real
AllenCahnPFFracture::computeQpResidual()
{
  return (-_gc_prop[_qp] * _l[_qp] * _grad_u[_qp] * _grad_test[_i][_qp] +
          2.0 * (1.0 - _u[_qp]) * _test[_i][_qp] * _H0_pos_old[_qp] -
          _gc_prop[_qp] / _l[_qp] * _u[_qp] * _test[_i][_qp]) /
         _gc_prop[_qp];
}

Real
AllenCahnPFFracture::computeQpJacobian()
{
  return (-_gc_prop[_qp] * _l[_qp] * _grad_phi[_j][_qp] * _grad_test[_i][_qp] -
          2.0 * _phi[_j][_qp] * _test[_i][_qp] * _H0_pos_old[_qp] -
          _gc_prop[_qp] / _l[_qp] * _phi[_j][_qp] * _test[_i][_qp]) /
         _gc_prop[_qp];
}

Real
AllenCahnPFFracture::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real u = _u[_qp];

  for (unsigned int c_comp = 0; c_comp < _ndisp; ++c_comp)
    if (jvar == _disp_var[c_comp])
    {
      if (u > 1.0)
        u = 1.0;
      const Real dxddFdc = 2.0 * (1.0 - u) * _test[_i][_qp];
      const Real d2Fdcdstrain_comp =
          (_d2Fdcdstrain[_qp].column(c_comp) + _d2Fdcdstrain[_qp].row(c_comp)) / 2.0 *
          _grad_phi[_j][_qp];
      return dxddFdc * d2Fdcdstrain_comp / _gc_prop[_qp];
    }

  return 0.0;
}
