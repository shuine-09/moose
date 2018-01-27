/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "XFEMTwoMaterialManagerConstraint.h"

#include "XFEMElemPairMaterialManager.h"
#include "MooseMesh.h"

template <>
InputParameters
validParams<XFEMTwoMaterialManagerConstraint>()
{
  InputParameters params = validParams<ElemElemConstraint>();
  params.addRequiredParam<UserObjectName>("manager", "XFEMElemPairMaterialManager object");
  params.addRequiredParam<UserObjectName>("manager_neighbor", "XFEMElemPairMaterialManager object");
  return params;
}

XFEMTwoMaterialManagerConstraint::XFEMTwoMaterialManagerConstraint(
    const InputParameters & parameters)
  : ElemElemConstraint(parameters),
    _manager(getUserObject<XFEMElemPairMaterialManager>("manager")),
    _manager_neighbor(getUserObject<XFEMElemPairMaterialManager>("manager_neighbor"))
{
}

void
XFEMTwoMaterialManagerConstraint::computeResidual()
{
  _manager.swapInProperties(std::min(_current_elem->unique_id(), _neighbor_elem->unique_id()));
  _manager_neighbor.swapInProperties(
      std::min(_current_elem->unique_id(), _neighbor_elem->unique_id()));
  ElemElemConstraint::computeResidual();
  _manager.swapOutProperties(std::min(_current_elem->unique_id(), _neighbor_elem->unique_id()));
  _manager_neighbor.swapOutProperties(
      std::min(_current_elem->unique_id(), _neighbor_elem->unique_id()));
}

void
XFEMTwoMaterialManagerConstraint::computeJacobian()
{
  _manager.swapInProperties(std::min(_current_elem->unique_id(), _neighbor_elem->unique_id()));
  _manager_neighbor.swapInProperties(
      std::min(_current_elem->unique_id(), _neighbor_elem->unique_id()));
  ElemElemConstraint::computeJacobian();
  _manager.swapOutProperties(std::min(_current_elem->unique_id(), _neighbor_elem->unique_id()));
  _manager_neighbor.swapOutProperties(
      std::min(_current_elem->unique_id(), _neighbor_elem->unique_id()));
}
