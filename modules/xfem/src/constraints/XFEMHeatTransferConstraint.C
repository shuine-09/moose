/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "XFEMHeatTransferConstraint.h"
#include "FEProblem.h"
#include "Assembly.h"

// libMesh includes
#include "libmesh/quadrature.h"

template<>
InputParameters validParams<XFEMHeatTransferConstraint>()
{
  InputParameters params = validParams<ElemElemConstraint>();
  params.set<Moose::MaterialDataType>("_material_data_type") = Moose::DIRAC_MATERIAL_DATA;
  return params;
}

XFEMHeatTransferConstraint::XFEMHeatTransferConstraint(const InputParameters & parameters) :
    ElemElemConstraint(parameters),
    MaterialPropertyInterface(this),
    _heat_flux(getMaterialProperty<Real>("heatflux"))
{
}

XFEMHeatTransferConstraint::~XFEMHeatTransferConstraint()
{
}

void
XFEMHeatTransferConstraint::reinitConstraintQuadrature(const ElementPairInfo & element_pair_info)
{
  _interface_normal = element_pair_info._elem1_normal;
  ElemElemConstraint::reinitConstraintQuadrature(element_pair_info);
}

Real
XFEMHeatTransferConstraint::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0;

  switch (type)
  {
    case Moose::Element:
      r += _heat_flux[_qp] * _test[_i][_qp];
      break;

    case Moose::Neighbor:
      r -= _heat_flux[_qp] * _test_neighbor[_i][_qp];
      break;
  }
  return r;
}

Real
XFEMHeatTransferConstraint::computeQpJacobian(Moose::DGJacobianType type)
{
  Real r = 0;

  return r;
}
