/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "XFEMResidualOpening.h"

// MOOSE includes
#include "Assembly.h"
#include "ElementPairInfo.h"
#include "FEProblem.h"

// // libMesh includes
// #include "libMesh/quadrature.h"

template <>
InputParameters
validParams<XFEMResidualOpening>()
{
  InputParameters params = validParams<XFEMMaterialManagerConstraint>();
  params.addParam<Real>("alpha", 1.0e5, "Penalty parameter.");
  params.addParam<Real>("roughness_fuel", 0.0, "Fuel roughness.");
  params.addParam<std::string>("base_name",
                               "Optional parameter that allows the user to define "
                               "multiple mechanics material systems on the same block");
  params.addCoupledVar("disp_x", "Coupled displacement in x");
  params.addCoupledVar("disp_y", "Coupled displacement in y");
  params.addRequiredParam<unsigned int>("component",
                                        "An integer corresponding to the direction "
                                        "the variable this kernel acts in. (0 for x, "
                                        "1 for y, 2 for z)");
  return params;
}

XFEMResidualOpening::XFEMResidualOpening(const InputParameters & parameters)
  : XFEMMaterialManagerConstraint(parameters),
    _disp_x(coupledValue("disp_x")),
    _disp_x_neighbor(coupledNeighborValue("disp_x")),
    _disp_y(coupledValue("disp_y")),
    _disp_y_neighbor(coupledNeighborValue("disp_y")),
    _component(getParam<unsigned int>("component")),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _alpha(getParam<Real>("alpha")),
    _roughness_fuel(getParam<Real>("roughness_fuel"))
{
}

void
XFEMResidualOpening::initialSetup()
{
  _max_normal_separation_old = getMaterialPropertyOld<Real>(_base_name + "max_normal_separation");
  _normal_separation = getMaterialPropertyOld<Real>(_base_name + "normal_separation");
}

XFEMResidualOpening::~XFEMResidualOpening() {}

void
XFEMResidualOpening::reinitConstraintQuadrature(const ElementPairInfo & element_pair_info)
{
  _interface_normal = element_pair_info._elem1_normal;
  ElemElemConstraint::reinitConstraintQuadrature(element_pair_info);
}

Real
XFEMResidualOpening::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0;
  const Real max_normal_separation = (*_max_normal_separation_old)[_qp];
  // const Real normal_separation = (*_normal_separation)[_qp];
  const Real normal_separation = (_disp_x_neighbor[_qp] - _disp_x[_qp]) * _interface_normal(0) +
                                 (_disp_y_neighbor[_qp] - _disp_y[_qp]) * _interface_normal(1);

  // Calculating Rotation matrix (from local to global coordinates):
  Real theta = std::atan2(_interface_normal(0), _interface_normal(1));
  Real R[2][2]; // Rotation matrix
  R[0][0] = std::cos(theta);
  R[1][0] = std::sin(theta);
  R[0][1] = -std::sin(theta);
  R[1][1] = std::cos(theta);

  // std::cout << "resdiual_opening = " << _roughness_fuel << std::endl;
  // std::cout << "max_normal_separation = " << max_normal_separation << std::endl;
  // std::cout << "normal_separation = " << normal_separation << std::endl;

  // Calculating Normal and tengential tractions on crack surface:
  Real t_n = 0.0;

  // if (max_normal_separation < _roughness_fuel)
  // {
  //   if (normal_separation < 0.0)
  //     t_n = _alpha * normal_separation;
  // }
  // else
  // {
  //   if (normal_separation < _roughness_fuel)
  //     t_n = _alpha * (normal_separation - _roughness_fuel);
  // }

  if (max_normal_separation < _roughness_fuel)
  {
    if (normal_separation < 0.5 * max_normal_separation)
      t_n = _alpha * (normal_separation - 0.5 * max_normal_separation);
  }
  else
  {
    if (normal_separation < 0.5 * _roughness_fuel)
      t_n = _alpha * (normal_separation - 0.5 * _roughness_fuel);
  }


 // if (normal_separation < max_normal_separation)
 //   t_n = _alpha * (normal_separation - max_normal_separation);

  // Rotating traction vector {t_n} to {t_x, t_y}:
  Real t_y = R[0][0] * t_n;
  Real t_x = R[1][0] * t_n;

  Real t_i = 0.0;
  if (_component == 0)
    t_i = t_x;
  else if (_component == 1)
    t_i = t_y;
  else
    mooseError("XFEMResidualOpening: 3D is not supported.");

  switch (type)
  {
    case Moose::Element:
      r -= t_i * _test[_i][_qp];
      break;

    case Moose::Neighbor:
      r += t_i * _test_neighbor[_i][_qp];
      break;
  }

  return r;
}

Real
XFEMResidualOpening::computeQpJacobian(Moose::DGJacobianType type)
{
  Real r = 0;
  const Real max_normal_separation = (*_max_normal_separation_old)[_qp];
  // const Real normal_separation = (*_normal_separation)[_qp];

  const Real normal_separation = (_disp_x_neighbor[_qp] - _disp_x[_qp]) * _interface_normal(0) +
                                 (_disp_y_neighbor[_qp] - _disp_y[_qp]) * _interface_normal(1);

  // Calculating Rotation matrix (from local to global coordinates):
  Real theta = std::atan2(_interface_normal(0), _interface_normal(1));
  Real R[2][2]; // Rotation matrix
  R[0][0] = std::cos(theta);
  R[1][0] = std::sin(theta);
  R[0][1] = -std::sin(theta);
  R[1][1] = std::cos(theta);

  Real factor = 0.0;

  // if (max_normal_separation < _roughness_fuel)
  // {
  //   if (normal_separation < 0.0)
  //     factor = _alpha;
  // }
  // else
  // {
  //   if (normal_separation < _roughness_fuel)
  //     factor = _alpha;
  // }

  if (max_normal_separation < _roughness_fuel)
  {
    if (normal_separation < 0.5 * max_normal_separation)
      factor = _alpha;
  }
  else
  {
    if (normal_separation < 0.5 * _roughness_fuel)
      factor = _alpha;
  }

//  if (normal_separation < max_normal_separation)
//     factor = _alpha;

  Real dseparation_du = 0.0;
  if (_component == 0)
    dseparation_du = R[1][0] * _interface_normal(0);
  else if (_component == 1)
    dseparation_du = R[0][0] * _interface_normal(1);
  else
    mooseError("XFEMResidualOpening: 3D is not supported.");

  switch (type)
  {
    case Moose::ElementElement:
      r -= factor * _test[_i][_qp] * _phi[_j][_qp] * dseparation_du;
      break;

    case Moose::ElementNeighbor:
      r += factor * _test[_i][_qp] * _phi_neighbor[_j][_qp] * dseparation_du;
      break;

    case Moose::NeighborElement:
      r += factor * _test_neighbor[_i][_qp] * _phi[_j][_qp] * dseparation_du;
      break;

    case Moose::NeighborNeighbor:
      r -= factor * _test_neighbor[_i][_qp] * _phi_neighbor[_j][_qp] * dseparation_du;
      break;
  }
  return r;
}
