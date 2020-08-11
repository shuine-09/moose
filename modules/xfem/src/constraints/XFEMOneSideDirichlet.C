//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "XFEMOneSideDirichlet.h"

// MOOSE includes
#include "Assembly.h"
#include "ElementPairInfo.h"
#include "FEProblem.h"
#include "GeometricCutUserObject.h"
#include "XFEM.h"

#include "libmesh/quadrature.h"

registerMooseObject("XFEMApp", XFEMOneSideDirichlet);

template <>
InputParameters
validParams<XFEMOneSideDirichlet>()
{
  InputParameters params = validParams<ElemElemConstraint>();
  params.addParam<Real>("alpha", 100, "Penalty parameter in penalty's formulation.");
  params.addRequiredParam<Real>("value",
                                "Dirichlet bc at specified side of the interface.");
  params.addRequiredParam<bool>("positive_side",
                                "Boolean specifying the side of the interface this constraint"
                                "applies to (positive side if true)");
  params.addParam<UserObjectName>(
      "geometric_cut_userobject",
      "Name of GeometricCutUserObject associated with this constraint.");
  return params;
}

XFEMOneSideDirichlet::XFEMOneSideDirichlet(const InputParameters & parameters)
  : ElemElemConstraint(parameters),
    _alpha(getParam<Real>("alpha")),
    _value(getParam<Real>("value")),
    _positive_side(getParam<bool>("positive_side"))
{
  _xfem = std::dynamic_pointer_cast<XFEM>(_fe_problem.getXFEM());
  if (_xfem == nullptr)
    mooseError("Problem casting to XFEM in XFEMOneSideDirichlet");

  const UserObject * uo =
      &(_fe_problem.getUserObjectBase(getParam<UserObjectName>("geometric_cut_userobject")));

  if (dynamic_cast<const GeometricCutUserObject *>(uo) == nullptr)
    mooseError("UserObject casting to GeometricCutUserObject in XFEMOneSideDirichlet");

  _interface_id = _xfem->getGeometricCutID(dynamic_cast<const GeometricCutUserObject *>(uo));
}

XFEMOneSideDirichlet::~XFEMOneSideDirichlet() {}

void
XFEMOneSideDirichlet::reinitConstraintQuadrature(const ElementPairInfo & element_pair_info)
{
  _interface_normal = element_pair_info._elem1_normal;
  ElemElemConstraint::reinitConstraintQuadrature(element_pair_info);
}

Real
XFEMOneSideDirichlet::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0;

  switch (type)
  {
    case Moose::Element:
      if (_positive_side)
      {
        if (_interface_normal(0) > 0.0)
          r += _alpha * (_u[_qp] - _value) * _test[_i][_qp];
      }
      else
      {
        if (_interface_normal(0) < 0.0)
          r += _alpha * (_u[_qp] - _value) * _test[_i][_qp];
      }
      break;

    case Moose::Neighbor:
      if (_positive_side)
      {
        if (_interface_normal(0) > 0.0)
          r += _alpha * (_u_neighbor[_qp] - _value) * _test_neighbor[_i][_qp];
      }
      else
      {
        if (_interface_normal(0) < 0.0)
          r += _alpha * (_u_neighbor[_qp] - _value) * _test_neighbor[_i][_qp];
      }
      break;
  }
  return r;
}

Real
XFEMOneSideDirichlet::computeQpJacobian(Moose::DGJacobianType type)
{
  Real r = 0;

  switch (type)
  {
    case Moose::ElementElement:
      r += _alpha * _phi[_j][_qp] * _test[_i][_qp];
      break;

    case Moose::NeighborNeighbor:
      r += _alpha * _phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp];
      break;
  }

  return r;
}
