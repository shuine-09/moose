//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "XFEMTwoSideDirichletNitsche.h"

// MOOSE includes
#include "Assembly.h"
#include "ElementPairInfo.h"
#include "FEProblem.h"
#include "GeometricCutUserObject.h"
#include "XFEM.h"

#include "libmesh/quadrature.h"

registerMooseObject("XFEMApp", XFEMTwoSideDirichletNitsche);

template <>
InputParameters
validParams<XFEMTwoSideDirichletNitsche>()
{
  InputParameters params = validParams<XFEMTwoSideDirichlet>();
  params.addRequiredParam<Real>("diffusivity_at_positive_level_set_side",
                                "Diffusivity at level set positive side.");
  params.addRequiredParam<Real>("diffusivity_at_negative_level_set_side",
                                "Diffusivity at level set negative side.");
  return params;
}

XFEMTwoSideDirichletNitsche::XFEMTwoSideDirichletNitsche(const InputParameters & parameters)
  : XFEMTwoSideDirichlet(parameters),
    _diffusivity_at_positive_level_set_side(
        getParam<Real>("diffusivity_at_positive_level_set_side")),
    _diffusivity_at_negative_level_set_side(
        getParam<Real>("diffusivity_at_negative_level_set_side"))
{
}

XFEMTwoSideDirichletNitsche::~XFEMTwoSideDirichletNitsche() {}

Real
XFEMTwoSideDirichletNitsche::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0;

  switch (type)
  {
    case Moose::Element:
      if (_interface_normal(0) > 0.0)
      {
        r += -_diffusivity_at_positive_level_set_side * _grad_u[_qp] * _interface_normal *
                 _test[_i][_qp] -
             (_value_at_positive_level_set_interface - _u[_qp]) *
                 _diffusivity_at_positive_level_set_side * _grad_test[_i][_qp] * _interface_normal;
        r += _alpha * (_u[_qp] - _value_at_positive_level_set_interface) * _test[_i][_qp];
      }
      else
      {
        r += -_diffusivity_at_negative_level_set_side * _grad_u[_qp] * _interface_normal *
                 _test[_i][_qp] -
             (_value_at_negative_level_set_interface - _u[_qp]) *
                 _diffusivity_at_negative_level_set_side * _grad_test[_i][_qp] * _interface_normal;
        r += _alpha * (_u[_qp] - _value_at_negative_level_set_interface) * _test[_i][_qp];
      }
      break;

    case Moose::Neighbor:
      if (_interface_normal(0) > 0.0)
      {
        r += -_diffusivity_at_negative_level_set_side * _grad_u_neighbor[_qp] * -_interface_normal *
                 _test_neighbor[_i][_qp] -
             (_value_at_negative_level_set_interface - _u_neighbor[_qp]) *
                 _diffusivity_at_negative_level_set_side * _grad_test_neighbor[_i][_qp] *
                 -_interface_normal;
        r += _alpha * (_u_neighbor[_qp] - _value_at_negative_level_set_interface) *
             _test_neighbor[_i][_qp];
      }
      else
      {
        r += -_diffusivity_at_positive_level_set_side * _grad_u_neighbor[_qp] * -_interface_normal *
                 _test_neighbor[_i][_qp] -
             (_value_at_positive_level_set_interface - _u_neighbor[_qp]) *
                 _diffusivity_at_positive_level_set_side * _grad_test_neighbor[_i][_qp] *
                 -_interface_normal;
        r += _alpha * (_u_neighbor[_qp] - _value_at_positive_level_set_interface) *
             _test_neighbor[_i][_qp];
      }
      break;
  }
  return r;
}

Real
XFEMTwoSideDirichletNitsche::computeQpJacobian(Moose::DGJacobianType type)
{
  Real r = 0;

  switch (type)
  {
    case Moose::ElementElement:
    {
      if (_interface_normal(0) > 0.0)
      {
        r += -_diffusivity_at_positive_level_set_side * _grad_phi[_j][_qp] * _interface_normal *
                 _test[_i][_qp] -
             (-_phi[_j][_qp]) * _diffusivity_at_positive_level_set_side * _grad_test[_i][_qp] *
                 _interface_normal;
        r += _alpha * _phi[_j][_qp] * _test[_i][_qp];
      }
      else
      {
        r += -_diffusivity_at_negative_level_set_side * _grad_phi[_j][_qp] * _interface_normal *
                 _test[_i][_qp] -
             (-_phi[_j][_qp]) * _diffusivity_at_negative_level_set_side * _grad_test[_i][_qp] *
                 _interface_normal;
        r += _alpha * _phi[_j][_qp] * _test[_i][_qp];
      }
    }
    break;

    case Moose::NeighborNeighbor:
    {
      if (_interface_normal(0) > 0.0)
      {
        r += -_diffusivity_at_negative_level_set_side * _grad_phi_neighbor[_j][_qp] *
                 -_interface_normal * _test_neighbor[_i][_qp] -
             (-_phi_neighbor[_j][_qp]) * _diffusivity_at_negative_level_set_side *
                 _grad_test_neighbor[_i][_qp] * -_interface_normal;
        r += _alpha * _phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp];
      }
      else
      {
        r += -_diffusivity_at_negative_level_set_side * _grad_phi_neighbor[_j][_qp] *
                 -_interface_normal * _test_neighbor[_i][_qp] -
             (-_phi_neighbor[_j][_qp]) * _diffusivity_at_negative_level_set_side *
                 _grad_test_neighbor[_i][_qp] * -_interface_normal;
        r += _alpha * _phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp];
      }
    }
    break;
  }

  return r;
}
