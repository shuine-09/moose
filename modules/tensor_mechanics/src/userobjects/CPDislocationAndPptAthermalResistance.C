/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "CPDislocationAndPptAthermalResistance.h"

registerMooseObject("TensorMechanicsApp", CPDislocationAndPptAthermalResistance);

template <>
InputParameters
validParams<CPDislocationAndPptAthermalResistance>()
{
  InputParameters params = validParams<CPDislocationBasedAthermalSlipResistance>();
  params.addParam<FunctionName>("precipitate_radius", 0.0, "Average radius of precipitate");
  params.addParam<Real>("precipitate_volume_fraction", 0.0, "Volume fraction of precipitate.");
  params.addParam<Real>("orowan_strength_factor", 1.0, "A prefactor for Orowan looping strength.");
  return params;
}

CPDislocationAndPptAthermalResistance::CPDislocationAndPptAthermalResistance(
    const InputParameters & parameters)
  : CPDislocationBasedAthermalSlipResistance(parameters),
    _precipitate_radius(getFunction("precipitate_radius")),
    _precipitate_volume_fraction(getParam<Real>("precipitate_volume_fraction")),
    _orowan_strength_factor(getParam<Real>("orowan_strength_factor"))
{
}

bool
CPDislocationAndPptAthermalResistance::calcSlipResistance(unsigned int qp,
                                                          std::vector<Real> & val) const
{
  const Real radius = 13e-3; //_precipitate_radius.value(_t, _q_point[qp]);

  if (CPDislocationBasedAthermalSlipResistance::calcSlipResistance(qp, val))
  {
    for (unsigned int i = 0; i < _variable_size; ++i)
    {
      Real disloc_resist = std::pow(val[i], 2.0);
      val[i] = 0.0;
      if (_precipitate_volume_fraction > 0.0)
      {
        Real spacing_precipitate =
            std::sqrt(8.0 / (libMesh::pi * 3.0 * _precipitate_volume_fraction)) * radius - radius;
        Real orowan_looping = _shear_mod * _b / spacing_precipitate * _orowan_strength_factor;
        val[i] += std::pow(orowan_looping, 2.0);
      }
      val[i] += disloc_resist;
      val[i] = std::sqrt(val[i]);
    }
  }
  else
    return false;

  return true;
}
