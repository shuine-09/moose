/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "XFEMNitscheDisplacementConstraint.h"

// MOOSE includes
#include "Assembly.h"
#include "ElementPairInfo.h"
#include "FEProblem.h"

// // libMesh includes
// #include "libMesh/quadrature.h"

template <>
InputParameters
validParams<XFEMNitscheDisplacementConstraint>()
{
  InputParameters params = validParams<XFEMTwoMaterialManagerConstraint>();
  params.addParam<std::string>("base_name",
                               "Optional parameter that allows the user to define "
                               "multiple mechanics material systems on the same block");
  params.addCoupledVar("disp_x", "Coupled displacement in x");
  params.addCoupledVar("disp_y", "Coupled displacement in y");
  params.addRequiredParam<unsigned int>("component",
                                        "An integer corresponding to the direction "
                                        "the variable this kernel acts in. (0 for x, "
                                        "1 for y, 2 for z)");
  params.addParam<Real>("alpha", 100, "Stablization parameter in Nitsche's formulation.");
  return params;
}

XFEMNitscheDisplacementConstraint::XFEMNitscheDisplacementConstraint(
    const InputParameters & parameters)
  : XFEMTwoMaterialManagerConstraint(parameters),
    _disp_x(coupledValue("disp_x")),
    _disp_x_neighbor(coupledNeighborValue("disp_x")),
    _disp_y(coupledValue("disp_y")),
    _disp_y_neighbor(coupledNeighborValue("disp_y")),
    _alpha(getParam<Real>("alpha")),
    _component(getParam<unsigned int>("component")),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : "")
{
}

void
XFEMNitscheDisplacementConstraint::initialSetup()
{
  _stress = getMaterialProperty<RankTwoTensor>(_base_name + "stress");
  _stress_neighbor = getNeighborMaterialProperty<RankTwoTensor>(_base_name + "stress");

  _elasticity_tensor = getMaterialProperty<RankFourTensor>(_base_name + "elasticity_tensor");
  _elasticity_tensor_neighbor =
      getNeighborMaterialProperty<RankFourTensor>(_base_name + "elasticity_tensor");
}

XFEMNitscheDisplacementConstraint::~XFEMNitscheDisplacementConstraint() {}

void
XFEMNitscheDisplacementConstraint::reinitConstraintQuadrature(
    const ElementPairInfo & element_pair_info)
{
  _interface_normal = element_pair_info._elem1_normal;
  ElemElemConstraint::reinitConstraintQuadrature(element_pair_info);
}

Real
XFEMNitscheDisplacementConstraint::computeQpResidual(Moose::DGResidualType type)
{
  RankTwoTensor stress;
  RankTwoTensor stress_neighbor;

  if (_current_elem->unique_id() ==
      std::min(_current_elem->unique_id(), _neighbor_elem->unique_id()))
  {
    stress = (*_stress)[_qp];
    stress_neighbor = (*_stress_neighbor)[_qp];
  }
  else
  {
    stress_neighbor = (*_stress)[_qp];
    stress = (*_stress_neighbor)[_qp];
  }

  Real u;
  Real u_neighbor;
  Real force = (stress * _interface_normal)(_component);
  Real force_neighbor = (stress_neighbor * _interface_normal)(_component);

  if (_component == 0)
  {
    u = _disp_x[_qp];
    u_neighbor = _disp_x_neighbor[_qp];
  }
  else if (_component == 1)
  {
    u = _disp_y[_qp];
    u_neighbor = _disp_y_neighbor[_qp];
  }

  RankTwoTensor grad_tensor(_grad_test[_i][_qp], _grad_test[_i][_qp], _grad_test[_i][_qp]);
  RankTwoTensor d_strain = (grad_tensor + grad_tensor.transpose()) / 2.0;
  RankTwoTensor d_stress = (*_elasticity_tensor)[_qp] * d_strain;

  RankTwoTensor grad_tensor_neighbor(
      _grad_test_neighbor[_i][_qp], _grad_test_neighbor[_i][_qp], _grad_test_neighbor[_i][_qp]);
  RankTwoTensor d_strain_neighbor = (grad_tensor_neighbor + grad_tensor_neighbor.transpose()) / 2.0;
  RankTwoTensor d_stress_neighbor = (*_elasticity_tensor_neighbor)[_qp] * d_strain_neighbor;

  Real r = 0;
  switch (type)
  {
    case Moose::Element:
      r -= (0.5 * force + 0.5 * force_neighbor) * _test[_i][_qp];
      r -= (u - u_neighbor) * 0.5 * (d_stress * _interface_normal)(_component);
      r += _alpha * (u - u_neighbor) * _test[_i][_qp];
      break;

    case Moose::Neighbor:
      r += (0.5 * force + 0.5 * force_neighbor) * _test_neighbor[_i][_qp];
      r -= (u - u_neighbor) * 0.5 * (d_stress_neighbor * _interface_normal)(_component);
      r -= _alpha * (u - u_neighbor) * _test_neighbor[_i][_qp];
      break;
  }
  return r;
}

Real
XFEMNitscheDisplacementConstraint::computeQpJacobian(Moose::DGJacobianType type)
{
  Real r = 0;

  RankTwoTensor grad_phi_tensor(_grad_phi[_j][_qp], _grad_phi[_j][_qp], _grad_phi[_j][_qp]);
  RankTwoTensor d_strain_phi = (grad_phi_tensor + grad_phi_tensor.transpose()) / 2.0;
  RankTwoTensor d_stress_phi = (*_elasticity_tensor)[_qp] * d_strain_phi;
  Real d_force_phi = (d_stress_phi * _interface_normal)(_component);

  RankTwoTensor grad_phi_tensor_neighbor(
      _grad_phi_neighbor[_j][_qp], _grad_phi_neighbor[_j][_qp], _grad_phi_neighbor[_j][_qp]);
  RankTwoTensor d_strain_phi_neighbor =
      (grad_phi_tensor_neighbor + grad_phi_tensor_neighbor.transpose()) / 2.0;
  RankTwoTensor d_stress_phi_neighbor = (*_elasticity_tensor_neighbor)[_qp] * d_strain_phi_neighbor;
  Real d_force_phi_neighbor = (d_stress_phi_neighbor * _interface_normal)(_component);

  RankTwoTensor grad_teset_tensor(_grad_test[_i][_qp], _grad_test[_i][_qp], _grad_test[_i][_qp]);
  RankTwoTensor d_strain_test = (grad_teset_tensor + grad_teset_tensor.transpose()) / 2.0;
  RankTwoTensor d_stress_test = (*_elasticity_tensor)[_qp] * d_strain_test;
  Real d_force_test = (d_stress_test * _interface_normal)(_component);

  RankTwoTensor grad_test_tensor_neighbor(
      _grad_test_neighbor[_i][_qp], _grad_test_neighbor[_i][_qp], _grad_test_neighbor[_i][_qp]);
  RankTwoTensor d_strain_test_neighbor =
      (grad_test_tensor_neighbor + grad_test_tensor_neighbor.transpose()) / 2.0;
  RankTwoTensor d_stress_test_neighbor =
      (*_elasticity_tensor_neighbor)[_qp] * d_strain_test_neighbor;
  Real d_force_test_neighbor = (d_stress_test_neighbor * _interface_normal)(_component);

  switch (type)
  {
    case Moose::ElementElement:
      r += -0.5 * d_force_phi * _test[_i][_qp] - _phi[_j][_qp] * 0.5 * d_force_test;
      r += _alpha * _phi[_j][_qp] * _test[_i][_qp];
      break;

    case Moose::ElementNeighbor:
      r += -0.5 * d_force_phi_neighbor * _test[_i][_qp] +
           _phi_neighbor[_j][_qp] * 0.5 * d_force_test;
      r -= _alpha * _phi_neighbor[_j][_qp] * _test[_i][_qp];
      break;

    case Moose::NeighborElement:
      r +=
          0.5 * d_force_phi * _test_neighbor[_i][_qp] - _phi[_j][_qp] * 0.5 * d_force_test_neighbor;
      r -= _alpha * _phi[_j][_qp] * _test_neighbor[_i][_qp];
      break;

    case Moose::NeighborNeighbor:
      r += 0.5 * d_force_phi_neighbor * _test_neighbor[_i][_qp] +
           _phi_neighbor[_j][_qp] * 0.5 * d_force_test_neighbor;
      r += _alpha * _phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp];
      break;
  }
  return r;
}
