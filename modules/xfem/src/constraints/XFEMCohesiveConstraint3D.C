/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "XFEMCohesiveConstraint3D.h"

// MOOSE includes
#include "Assembly.h"
#include "ElementPairInfo.h"
#include "FEProblem.h"

// // libMesh includes
// #include "libMesh/quadrature.h"

template <>
InputParameters
validParams<XFEMCohesiveConstraint3D>()
{
  InputParameters params = validParams<XFEMMaterialManagerConstraint>();
  params.addParam<Real>("stiffness", 9.6628, "Initial stiffness.");
  params.addParam<Real>("max_traction", 11.75, "Max traction.");
  params.addParam<std::string>("base_name",
                               "Optional parameter that allows the user to define "
                               "multiple mechanics material systems on the same block");
  params.addParam<Real>("Gc", 82.2478, "Strain energy release rate.");
  params.addCoupledVar("disp_x", "Coupled displacement in x");
  params.addCoupledVar("disp_y", "Coupled displacement in y");
  params.addCoupledVar("disp_z", "Coupled displacement in y");
  params.addRequiredParam<unsigned int>("component",
                                        "An integer corresponding to the direction "
                                        "the variable this kernel acts in. (0 for x, "
                                        "1 for y, 2 for z)");
  return params;
}

XFEMCohesiveConstraint3D::XFEMCohesiveConstraint3D(const InputParameters & parameters)
  : XFEMMaterialManagerConstraint(parameters),
    _stiffness(getParam<Real>("stiffness")),
    _max_traction(getParam<Real>("max_traction")),
    _Gc(getParam<Real>("Gc")),
    _disp_x(coupledValue("disp_x")),
    _disp_x_neighbor(coupledNeighborValue("disp_x")),
    _disp_y(coupledValue("disp_y")),
    _disp_y_neighbor(coupledNeighborValue("disp_y")),
    _disp_z(coupledValue("disp_z")),
    _disp_z_neighbor(coupledNeighborValue("disp_z")),
    _component(getParam<unsigned int>("component")),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : "")
{
}

void
XFEMCohesiveConstraint3D::initialSetup()
{
  _max_normal_separation_old = getMaterialPropertyOld<Real>(_base_name + "max_normal_separation");
}

XFEMCohesiveConstraint3D::~XFEMCohesiveConstraint3D() {}

void
XFEMCohesiveConstraint3D::reinitConstraintQuadrature(const ElementPairInfo & element_pair_info)
{
  _interface_normal = element_pair_info._elem1_normal;
  ElemElemConstraint::reinitConstraintQuadrature(element_pair_info);
}

Real
XFEMCohesiveConstraint3D::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0;
  Real delta_0 = _max_traction / _stiffness; // separation at max traction
  Real delta_f = 2. * _Gc / _max_traction;   // separation at failure for bilinear form
  const Real max_normal_separation = (*_max_normal_separation_old)[_qp];

  // Calculating normal separation:
  Real delta_m_n = (_disp_x_neighbor[_qp] - _disp_x[_qp]) * _interface_normal(0) +
                   (_disp_y_neighbor[_qp] - _disp_y[_qp]) * _interface_normal(1) +
                   (_disp_z_neighbor[_qp] - _disp_z[_qp]) * _interface_normal(2);

  Real stiffness_deg = _stiffness; // Initializing stiffness
  Real alpha = 1.0e6;

  // Calculating Normal and tengential tractions on crack surface:
  Real t_n = 0.0;

  // std::cout << "separation = " << delta_m_n << ", max_normal_separation = " <<
  // max_normal_separation
  //           << ", delta_0 = " << delta_0 << ", delta_f = " << delta_f << std::endl;

  // std::cout << "delta_0 = " << delta_0 << ", delta_f = " << delta_f
  //           << ", max_normal_separation = " << max_normal_separation << std::endl;

  if (max_normal_separation < delta_0)
  {
    // std::cout << "max_normal_separation < delta_0"
    //           << ", delta_m_n = " << delta_m_n << std::endl;
    // std::cout << "max_normal_separation = " << max_normal_separation << ", delta_0 = " << delta_0
    //           << std::endl;
    if (delta_m_n < 0.0)
      t_n = alpha * delta_m_n;
    else if (delta_m_n >= 0.0 && delta_m_n <= delta_0)
      t_n = _stiffness * delta_m_n;
    else if (delta_m_n > delta_0 && delta_m_n < delta_f)
      t_n = (_max_traction / (delta_f - delta_0)) * (delta_f - delta_m_n);
    else if (delta_m_n >= delta_f)
      t_n = 0.0;
  }
  else if (max_normal_separation >= delta_0 && max_normal_separation < delta_f)
  {
    // std::cout << "max_normal_separation >= delta_0 && max_normal_separation < delta_f" <<
    // std::endl;
    if (delta_m_n < 0.0)
      t_n = alpha * (delta_m_n - 0.0);
    else if (delta_m_n >= 0.0 && delta_m_n <= max_normal_separation)
    {
      stiffness_deg = (_max_traction / (delta_f - delta_0)) * (delta_f - max_normal_separation) /
                      max_normal_separation;
      t_n = stiffness_deg * delta_m_n;
    }
    else if (delta_m_n > max_normal_separation && delta_m_n < delta_f)
    {
      t_n = (_max_traction / (delta_f - delta_0)) * (delta_f - delta_m_n);
    }
    else if (delta_m_n >= delta_f)
      t_n = 0.0;
  }
  else if (max_normal_separation >= delta_f)
  {
    // std::cout << "max_normal_separation >= delta_f" << std::endl;
    t_n = 0.0;
  }

  // std::cout << "max_normal_separation = " << max_normal_separation << ", delta_m_n = " <<
  // delta_m_n
  //           << ", tn = " << t_n << std::endl;

  // Rotating traction vector {t_n} to {t_x, t_y, t_z}:
  Real t_x = _interface_normal(0) * t_n;
  Real t_y = _interface_normal(1) * t_n;
  Real t_z = _interface_normal(2) * t_n;

  // std::cout << ">>>>>>>>>>>>>>>>>>>>>>>> " << std::endl;
  // std::cout << "t_x = " << t_x << ", t_y = " << t_y << std::endl;
  // std::cout << "R1 = " << R[0][0] << ", " << R[0][1] << std::endl;
  // std::cout << "R2 = " << R[1][0] << ", " << R[1][1] << std::endl;

  Real t_i = 0.0;
  if (_component == 0)
    t_i = t_x;
  else if (_component == 1)
    t_i = t_y;
  else if (_component == 2)
    t_i = t_z;
  else
    mooseError("XFEMCohesiveConstraint3D: 3D is not supported.");

  switch (type)
  {
    case Moose::Element:
      r -= t_i * _test[_i][_qp];
      break;

    case Moose::Neighbor:
      r += t_i * _test_neighbor[_i][_qp];
      break;
  }
  // std::cout << "t_i = " << t_i << std::endl;
  // std::cout << "normal = " << _interface_normal << std::endl;
  // std::cout << "delta_m_n = " << delta_m_n << ", r[" << _component << "] = " << r << std::endl;
  // std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
  // const_cast<ComputeCohesiveTraction &>(_comp_trac).setSeparation(delta_m_n);
  // _traction[_qp] = r;
  return r;
}

Real
XFEMCohesiveConstraint3D::computeQpJacobian(Moose::DGJacobianType type)
{
  Real r = 0;

  Real delta_0 = _max_traction / _stiffness; // separation at max traction
  Real delta_f = 2. * _Gc / _max_traction;   // separation at failure for bilinear form
  const Real max_normal_separation = (*_max_normal_separation_old)[_qp];

  // Calculating normal separation:
  Real delta_m_n = (_disp_x_neighbor[_qp] - _disp_x[_qp]) * _interface_normal(0) +
                   (_disp_y_neighbor[_qp] - _disp_y[_qp]) * _interface_normal(1) +
                   (_disp_z_neighbor[_qp] - _disp_z[_qp]) * _interface_normal(2);

  // Point interface_tangent(_interface_normal(1), -_interface_normal(0), 0);
  //
  // // Calculating tangential separation:
  // Real delta_m_t = (_disp_x_neighbor[_qp] - _disp_x[_qp]) * interface_tangent(0) +
  //                  (_disp_y_neighbor[_qp] - _disp_y[_qp]) * interface_tangent(1);

  Real factor = _stiffness; // jacobian factor
  Real alpha = 1.0e6;

  if (max_normal_separation < delta_0)
  {
    if (delta_m_n < 0.0)
      factor = alpha;
    else if (delta_m_n >= 0.0 && delta_m_n <= delta_0)
      factor = _stiffness;
    else if (delta_m_n > delta_0 && delta_m_n < delta_f)
      factor = -(_max_traction / (delta_f - delta_0));
    else if (delta_m_n >= delta_f)
      factor = 0.0;
  }
  else if (max_normal_separation >= delta_0 && max_normal_separation < delta_f)
  {
    if (delta_m_n < 0.0)
      factor = alpha;
    else if (delta_m_n >= 0.0 && delta_m_n <= max_normal_separation)
      factor = (_max_traction / (delta_f - delta_0)) * (delta_f - max_normal_separation) /
               max_normal_separation;
    else if (delta_m_n > max_normal_separation && delta_m_n < delta_f)
      factor = -(_max_traction / (delta_f - delta_0));
    else if (delta_m_n >= delta_f)
      factor = 0.0;
  }
  else if (max_normal_separation >= delta_f)
    factor = 0.0;

  Real dseparation_du = 0.0;
  if (_component == 0)
    dseparation_du = _interface_normal(0) * _interface_normal(0);
  else if (_component == 1)
    dseparation_du = _interface_normal(1) * _interface_normal(1);
  else if (_component == 2)
    dseparation_du = _interface_normal(2) * _interface_normal(2);
  else
    mooseError("XFEMCohesiveConstraint3D: 3D is not supported.");

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
