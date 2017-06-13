/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "CPDislocationBasedDensityOrowanLoopingRateComp.h"

registerMooseObject("TensorMechanicsApp", CPDislocationBasedDensityOrowanLoopingRateComp);

template <>
InputParameters
validParams<CPDislocationBasedDensityOrowanLoopingRateComp>()
{
  InputParameters params = validParams<CrystalPlasticityStateVarRateComponent>();
  params.addParam<std::string>("uo_glide_slip_rate_name",
                               "Name of glide slip rate: Same as state "
                               "variable user object specified in input file.");
  params.addRequiredParam<Real>("r_cl", "Radius corresponding to the loss of coherence");
  params.addRequiredParam<Real>(
      "r_trans", "Radius corrsponding the transition from shearable to non-shearable");
  params.addRequiredParam<Real>("factor", "Radius corresponding to the loss of coherence");
  params.addParam<Real>("precipitate_volume_fraction", 0.0, "Volume fraction of precipitate.");
  params.addParam<Real>("precipitate_radius", 0.0, "Average radius of precipitate");

  params.addClassDescription("Dislocation based constitutive model userobject "
                             "class for mobile dislocation density rate "
                             "component for orowan looping mechanism. Override the "
                             "virtual functions in your class");
  return params;
}

CPDislocationBasedDensityOrowanLoopingRateComp::CPDislocationBasedDensityOrowanLoopingRateComp(
    const InputParameters & parameters)
  : CrystalPlasticityStateVarRateComponent(parameters),
    _mat_prop_glide_slip_rate(getMaterialProperty<std::vector<Real>>(
        parameters.get<std::string>("uo_glide_slip_rate_name"))),
    _r_cl(getParam<Real>("r_cl")),
    _r_trans(getParam<Real>("r_trans")),
    _factor(getParam<Real>("factor")),
    _precipitate_radius(getParam<Real>("precipitate_radius")),
    _precipitate_volume_fraction(getParam<Real>("precipitate_volume_fraction"))
{
}

bool
CPDislocationBasedDensityOrowanLoopingRateComp::calcStateVariableEvolutionRateComponent(
    unsigned int qp, std::vector<Real> & val) const
{
  val.assign(_variable_size, 0.0);

  Real efficiency = 0.0;

  if (_precipitate_radius > _r_cl)
    efficiency = 1.0;
  else if (_precipitate_radius > _r_trans)
    efficiency = (_precipitate_radius - _r_trans) / (_r_cl - _r_trans);
  else
    efficiency = 0.0;

  for (unsigned int i = 0; i < _variable_size; ++i)
  {
    val[i] = efficiency * _factor * std::sqrt(_precipitate_volume_fraction) / _precipitate_radius *
             std::abs(_mat_prop_glide_slip_rate[qp][i]);
  }

  return true;
}
