//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputePFFractureStressBase.h"

defineLegacyParams(ComputePFFractureStressBase);

InputParameters
ComputePFFractureStressBase::validParams()
{
  InputParameters params = ComputeStressBase::validParams();
  params.addRequiredCoupledVar("c", "Name of damage variable");
  params.addParam<bool>(
      "use_current_history_variable", false, "Use the current value of the history variable.");
  params.addParam<bool>("use_snes_vi_solver",
                        false,
                        "Use PETSc's SNES variational inequalities solver to enforce damage "
                        "irreversibility condition and restrict damage between 0 and 1.");
  params.addParam<MaterialPropertyName>("barrier_energy",
                                        "Name of material property for fracture energy barrier.");
  params.addParam<MaterialPropertyName>(
      "E_name", "elastic_energy", "Name of material property for elastic energy");
  params.addParam<MaterialPropertyName>(
      "D_name", "degradation", "Name of material property for energetic degradation function.");
  params.addParam<MaterialPropertyName>(
      "I_name", "indicator", "Name of material property for damage indicator function.");
  params.addParam<MaterialPropertyName>(
      "F_name",
      "local_fracture_energy",
      "Name of material property for local fracture energy function.");
  return params;
}

ComputePFFractureStressBase::ComputePFFractureStressBase(const InputParameters & parameters)
  : ComputeStressBase(parameters),
    _elasticity_tensor_name(_base_name + "elasticity_tensor"),
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_elasticity_tensor_name)),
    _c(coupledValue("c")),
    _l(getMaterialProperty<Real>("l")),
    _gc(getMaterialProperty<Real>("gc_prop")),
    _pressure(getDefaultMaterialProperty<Real>("fracture_pressure")),
    _use_current_hist(getParam<bool>("use_current_history_variable")),
    _use_snes_vi_solver(getParam<bool>("use_snes_vi_solver")),
    _H(declareProperty<Real>("hist")),
    _H_old(getMaterialPropertyOld<Real>("hist")),
    _barrier(getDefaultMaterialProperty<Real>("barrier_energy")),
    _E(declareProperty<Real>(getParam<MaterialPropertyName>("E_name"))),
    _dEdc(declarePropertyDerivative<Real>(getParam<MaterialPropertyName>("E_name"),
                                          getVar("c", 0)->name())),
    _d2Ed2c(declarePropertyDerivative<Real>(
        getParam<MaterialPropertyName>("E_name"), getVar("c", 0)->name(), getVar("c", 0)->name())),
    _dstress_dc(
        declarePropertyDerivative<RankTwoTensor>(_base_name + "stress", getVar("c", 0)->name())),
    _d2Fdcdstrain(declareProperty<RankTwoTensor>("d2Fdcdstrain")),
    _D(getMaterialProperty<Real>("D_name")),
    _dDdc(getMaterialPropertyDerivative<Real>("D_name", getVar("c", 0)->name())),
    _d2Dd2c(getMaterialPropertyDerivative<Real>(
        "D_name", getVar("c", 0)->name(), getVar("c", 0)->name())),
    _I(getDefaultMaterialProperty<Real>("I_name")),
    _dIdc(getMaterialPropertyDerivative<Real>("I_name", getVar("c", 0)->name())),
    _d2Id2c(getMaterialPropertyDerivative<Real>(
        "I_name", getVar("c", 0)->name(), getVar("c", 0)->name()))
{
}

void
ComputePFFractureStressBase::initQpStatefulProperties()
{
  _H[_qp] = 0.0;

  // Real d;
  // Real x = _q_point[_qp](0);
  // Real y = _q_point[_qp](1);
  //
  // Real x_center = 1.8;
  // Real y_center = 2;
  //
  // Real x_center2 = 2.2;
  // Real y_center2 = 2;
  //
  // // if (x < x_center)
  // //   d = std::abs(y - y_center);
  // // else
  // //   d = std::sqrt((x - x_center) * (x - x_center) + (y - y_center) * (y - y_center));
  //
  // if (x < x_center2 && x > x_center)
  //   d = std::abs(y - y_center);
  // else if (x < x_center)
  //   d = std::sqrt((x - x_center) * (x - x_center) + (y - y_center) * (y - y_center));
  // else
  //   d = std::sqrt((x - x_center2) * (x - x_center2) + (y - y_center2) * (y - y_center2));
  //
  // if (d <= _l[_qp])
  //   _H[_qp] = 100 * (_gc[_qp] / 2.0 / _l[_qp] * (1 - 2.0 * d / _l[_qp]));
  // else
  //   _H[_qp] = 0;
}
