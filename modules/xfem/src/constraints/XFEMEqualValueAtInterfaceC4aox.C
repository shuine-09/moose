//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "XFEMEqualValueAtInterfaceC4aox.h"
#include "FEProblem.h"
#include "GeometricCutUserObject.h"
#include "XFEM.h"

registerMooseObject("XFEMApp", XFEMEqualValueAtInterfaceC4aox);

InputParameters
XFEMEqualValueAtInterfaceC4aox::validParams()
{
  InputParameters params = ElemElemConstraint::validParams();
  params.addRequiredParam<Real>("alpha", "Penalty parameter in penalty formulation.");
  params.addRequiredParam<Real>("temperature", "Temperature [K] at the alpha/beta the interface");
  params.addParam<bool>("offset",true, "Specifies if an offset is needed or not. Due to the use "
                                "of a weak discontinuity variable, an offset on the oxide/alpha interface"
                                "value is needed when an alpha/beta interface is also used.");
  params.addParam<UserObjectName>(
      "geometric_cut_userobject",
      "Name of GeometricCutUserObject associated with this constraint.");
  params.addClassDescription("Enforce that the solution have the same value on opposing sides of "
                             "an XFEM interface. Here this is the oxide/alpha interface value is given"
                             "as modelled in the C4 model. The value is given by Zr phase diagram.");
  return params;
}

XFEMEqualValueAtInterfaceC4aox::XFEMEqualValueAtInterfaceC4aox(const InputParameters & parameters)
  : ElemElemConstraint(parameters), _alpha(getParam<Real>("alpha")), _temperature(getParam<Real>("temperature")),
    _offset(getParam<bool>("offset"))
{
  _xfem = std::dynamic_pointer_cast<XFEM>(_fe_problem.getXFEM());
  if (_xfem == nullptr)
    mooseError("Problem casting to XFEM in XFEMEqualValueAtInterfaceC4aox");

  const UserObject * uo =
      &(_fe_problem.getUserObjectBase(getParam<UserObjectName>("geometric_cut_userobject")));

  if (dynamic_cast<const GeometricCutUserObject *>(uo) == nullptr)
    mooseError("UserObject casting to GeometricCutUserObject in XFEMEqualValueAtInterfaceC4aox");

  _interface_id = _xfem->getGeometricCutID(dynamic_cast<const GeometricCutUserObject *>(uo));

  Real x_o_a_ox = 0.2978;
  if (_temperature > 473.15 && _temperature < 1478.15)
  {
    x_o_a_ox = (28.6 + exp(-6748/_temperature + 4.748)) * 1e-2;
  }
  else if (_temperature > 1478.15 && _temperature < 1798.15)
  {
    x_o_a_ox = (28.6 + exp(-6301/_temperature + 4.460)) * 1e-2;
  }
  else if (_temperature > 1798.15 && _temperature < 2338.15)
  {
    x_o_a_ox = (28.6 + exp(-7012/_temperature + 8.434 - 3.521e-3 * _temperature )) * 1e-2;
  }
  else
  {
    x_o_a_ox = 28.6 *1e-2;
  }

  _value = x_o_a_ox / (1 - x_o_a_ox);

  if (_offset)
  {
    const Real x_o_b_a = (9.59e-3 * (_temperature - 1136) + 4.72e-6 * pow(_temperature - 1136,2) - 4.35e-9 * pow(_temperature - 1136,3)) * 1e-2;
    const Real x_o_a_b = (45.86e-3 * (_temperature - 1136) - 44.77e-6 * pow(_temperature - 1136,2) + 17.40e-9 * pow(_temperature - 1136,3)) * 1e-2; //the original one, not the weak equivalent
    const Real c_o_a_b = x_o_a_b / (1 - x_o_a_b);
    const Real c_o_b_a = x_o_b_a / (1 - x_o_b_a);
    _value = _value - (c_o_a_b - c_o_b_a);
  }


}

XFEMEqualValueAtInterfaceC4aox::~XFEMEqualValueAtInterfaceC4aox() {}

void
XFEMEqualValueAtInterfaceC4aox::reinitConstraintQuadrature(const ElementPairInfo & element_pair_info)
{
  ElemElemConstraint::reinitConstraintQuadrature(element_pair_info);
}

Real
XFEMEqualValueAtInterfaceC4aox::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0;

  switch (type)
  {
    case Moose::Element:
      r += _alpha * (_u[_qp] - _value) * _test[_i][_qp];
      break;

    case Moose::Neighbor:
      r += _alpha * (_u_neighbor[_qp] - _value) * _test_neighbor[_i][_qp];
      break;
  }
  return r;
}

Real
XFEMEqualValueAtInterfaceC4aox::computeQpJacobian(Moose::DGJacobianType type)
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

    default:
      break;
  }

  return r;
}
