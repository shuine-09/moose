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

#include "VesicleVolumeAreaPenalty.h"

registerMooseObject("MooseApp", VesicleVolumeAreaPenalty);

template <>
InputParameters
validParams<VesicleVolumeAreaPenalty>()
{
  InputParameters params = validParams<Kernel>();
  params.addParam<Real>("alpha_v", 1000, "The penalty parameter of vesicle volume.");
  params.addParam<Real>("alpha_a", 1000, "The penalty parameter of vesicle area.");
  params.addParam<Real>("epsilon", 0.01, "The interfacial penalty parameter.");
  params.addParam<bool>("rz", false, "The RZ coordindate system.");
  params.addRequiredParam<PostprocessorName>("vesicle_volume",
                                             "Name of vesicle volume user object.");
  params.addRequiredParam<PostprocessorName>("vesicle_area", "Name of vesicle area user object.");
  params.addParam<bool>("use_prescribed_volume", false, "Use prescribed volume.");
  params.addParam<bool>("use_prescribed_area", false, "Use prescribed area.");
  params.addParam<Real>("prescribed_volume", 0.0, "Prescribed volume.");
  params.addParam<Real>("prescribed_area", 0.0, "Prescribed area.");
  params.addParam<bool>("use_nonlocal_constraint", false, "Use non-local constraint.");
  return params;
}

VesicleVolumeAreaPenalty::VesicleVolumeAreaPenalty(const InputParameters & parameters)
  : Kernel(parameters),
    _alpha_v(getParam<Real>("alpha_v")),
    _alpha_a(getParam<Real>("alpha_a")),
    _epsilon(getParam<Real>("epsilon")),
    _rz(getParam<bool>("rz")),
    _second_phi(secondPhi()),
    _second_test(secondTest()),
    _second_u(second()),
    _vesicle_area(getPostprocessorValue("vesicle_area")),
    _vesicle_volume(getPostprocessorValue("vesicle_volume")),
    _use_prescribed_volume(getParam<bool>("use_prescribed_volume")),
    _use_prescribed_area(getParam<bool>("use_prescribed_area")),
    _prescribed_volume(getParam<Real>("prescribed_volume")),
    _prescribed_area(getParam<Real>("prescribed_area")),
    _use_nonlocal_constraint(getParam<bool>("use_nonlocal_constraint"))
{
  _alpha_v0 = _alpha_v;
  _alpha_a0 = _alpha_a;
}

VesicleVolumeAreaPenalty::~VesicleVolumeAreaPenalty() {}

void
VesicleVolumeAreaPenalty::timestepSetup()
{
  if (_t_step <= 5)
  {
    _volume_0 = _vesicle_volume;
    _area_0 = _vesicle_area;
    _alpha_v = 0.0;
    _alpha_a = 0.0;
  }
  else
  {
    _alpha_a = _alpha_a0;
    _alpha_v = _alpha_v0;
  }

  if (_use_prescribed_volume)
  {
    _volume_0 = _prescribed_volume;
    if (_t_step <= 500)
      _alpha_v = _alpha_v0 / 500 * _t_step;
    else
      _alpha_v = _alpha_v0;
  }

  if (_use_prescribed_area)
  {
    _area_0 = _prescribed_area;
    if (_t_step <= 500)
      _alpha_a = _alpha_a0 / 500 * _t_step;
    else
      _alpha_a = _alpha_a0;
  }
}

void
VesicleVolumeAreaPenalty::jacobianSetup()
{
  _volume = _vesicle_volume;
  _area = _vesicle_area;
}

void
VesicleVolumeAreaPenalty::residualSetup()
{
  // std::cout << "volume = " << _volume << ", volume_0 = " << _volume_0
  //           << ", diff = " << _volume - _volume_0 << std::endl;
  // std::cout << "area = " << _area << ", area_0 = " << _area_0 << ", diff= " << _area - _area_0
  //           << std::endl;
  _volume = _vesicle_volume;
  _area = _vesicle_area;
}

Real
VesicleVolumeAreaPenalty::computeQpResidual()
{
  Real r = 0;

  Real rz_coord = _q_point[_qp](0);

  Real lap_u = _second_u[_qp].tr();

  if (_rz)
    lap_u += _grad_u[_qp](0) / rz_coord;

  if (!_use_nonlocal_constraint)
  {
    r += -_alpha_v * _test[_i][_qp] * (_volume - _volume_0);

    r += _alpha_a * _test[_i][_qp] * (_area - _area_0) *
         (-3.0 / 2.0 / std::sqrt(2.0) * _epsilon *
          (lap_u - 1.0 / pow(_epsilon, 2.0) * (pow(_u[_qp], 2.0) - 1.0) * _u[_qp]));
  }
  else
  {
    r += -_alpha_v * _test[_i][_qp] * (-_volume_0);

    r += _alpha_a * _test[_i][_qp] * (-_area_0) * (-3.0 / std::sqrt(2.0) * _epsilon * lap_u);
  }

  return r;
}

Real
VesicleVolumeAreaPenalty::computeQpJacobian()
{
  Real r = 0;

  Real rz_coord = _q_point[_qp](0);

  Real lap_phi_j = _second_phi[_j][_qp].tr();

  if (_rz)
    lap_phi_j += _grad_phi[_j][_qp](0) / rz_coord;

  r += _alpha_a * _test[_i][_qp] * (_area - _area_0) *
       (-3.0 / 2.0 / std::sqrt(2.0) * _epsilon *
        (lap_phi_j -
         1.0 / pow(_epsilon, 2.0) * (3.0 * pow(_u[_qp], 2.0) * _phi[_j][_qp] - _phi[_j][_qp])));

  return r;
}
