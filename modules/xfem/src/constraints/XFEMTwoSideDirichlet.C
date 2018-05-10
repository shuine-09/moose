//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "XFEMTwoSideDirichlet.h"

// MOOSE includes
#include "Assembly.h"
#include "ElementPairInfo.h"
#include "FEProblem.h"
#include "GeometricCutUserObject.h"
#include "XFEM.h"

#include "libmesh/quadrature.h"

registerMooseObject("XFEMApp", XFEMTwoSideDirichlet);

template <>
InputParameters
validParams<XFEMTwoSideDirichlet>()
{
  InputParameters params = validParams<ElemElemConstraint>();
  params.addParam<Real>("alpha", 100, "Stablization parameter in Nitsche's formulation.");
  params.addRequiredParam<Real>("level_set_positive_value",
                                "Dirichlet bc for level set positive region.");
  params.addRequiredParam<Real>("level_set_negative_value",
                                "Dirichlet bc for level set negative region.");
  params.addParam<UserObjectName>(
      "geometric_cut_userobject",
      "Name of GeometricCutUserObject associated with this constraint.");
  return params;
}

XFEMTwoSideDirichlet::XFEMTwoSideDirichlet(const InputParameters & parameters)
  : ElemElemConstraint(parameters),
    _alpha(getParam<Real>("alpha")),
    _levelset_positive_value(getParam<Real>("level_set_positive_value")),
    _levelset_negative_value(getParam<Real>("level_set_negative_value"))
{
  _xfem = std::dynamic_pointer_cast<XFEM>(_fe_problem.getXFEM());
  if (_xfem == nullptr)
    mooseError("Problem casting to XFEM in XFEMTwoSideDirichlet");

  const UserObject * uo =
      &(_fe_problem.getUserObjectBase(getParam<UserObjectName>("geometric_cut_userobject")));

  if (dynamic_cast<const GeometricCutUserObject *>(uo) == nullptr)
    mooseError("UserObject casting to GeometricCutUserObject in XFEMTwoSideDirichlet");

  _interface_id = _xfem->getGeometricCutID(dynamic_cast<const GeometricCutUserObject *>(uo));
}

XFEMTwoSideDirichlet::~XFEMTwoSideDirichlet() {}

void
XFEMTwoSideDirichlet::reinitConstraintQuadrature(const ElementPairInfo & element_pair_info)
{
  _interface_normal = element_pair_info._elem1_normal;
  ElemElemConstraint::reinitConstraintQuadrature(element_pair_info);
}

Real
XFEMTwoSideDirichlet::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0;

  switch (type)
  {
    case Moose::Element:
      if (_interface_normal(0) > 0.0)
        r += _alpha * (_u[_qp] - _levelset_positive_value) * _test[_i][_qp];
      else
        r += _alpha * (_u[_qp] - _levelset_negative_value) * _test[_i][_qp];
      break;

    case Moose::Neighbor:
      if (_interface_normal(0) > 0.0)
        r += _alpha * (_u_neighbor[_qp] - _levelset_negative_value) * _test_neighbor[_i][_qp];
      else
        r += _alpha * (_u_neighbor[_qp] - _levelset_positive_value) * _test_neighbor[_i][_qp];
      break;
  }
  return r;
}

Real
XFEMTwoSideDirichlet::computeQpJacobian(Moose::DGJacobianType type)
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
