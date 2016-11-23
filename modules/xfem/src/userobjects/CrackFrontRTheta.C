/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "CrackFrontRTheta.h"
#include "libmesh/fe_interface.h"
#include "DisplacedProblem.h"

template<>
InputParameters validParams<CrackFrontRTheta>()
{
  InputParameters params = validParams<DiscreteElementUserObject>();
  return params;
}

CrackFrontRTheta::CrackFrontRTheta(const InputParameters & parameters):
    DiscreteElementUserObject(parameters)
{
}

RealVectorValue
CrackFrontRTheta::rotateToCrackFrontCoords(const RealVectorValue vector) const
{
  ColumnMajorMatrix vec3x1;
  vec3x1 = _rot_mat * vector;
  RealVectorValue vec;
  vec(0) = vec3x1(0,0);
  vec(1) = vec3x1(1,0);
  vec(2) = vec3x1(2,0);
  return vec;
}

RealVectorValue
CrackFrontRTheta::rotateToCrackFrontCoords(const Point point) const
{
  RealVectorValue vector(point(0), point(1), point(2));
  ColumnMajorMatrix vec3x1;
  vec3x1 = _rot_mat * vector;
  RealVectorValue vec;
  vec(0) = vec3x1(0,0);
  vec(1) = vec3x1(1,0);
  vec(2) = vec3x1(2,0);
  return vec;
}


ColumnMajorMatrix
CrackFrontRTheta::rotateToCrackFrontCoords(const SymmTensor tensor) const
{
  ColumnMajorMatrix tensor_CMM;
  tensor_CMM(0,0) = tensor.xx();
  tensor_CMM(0,1) = tensor.xy();
  tensor_CMM(0,2) = tensor.xz();
  tensor_CMM(1,0) = tensor.xy();
  tensor_CMM(1,1) = tensor.yy();
  tensor_CMM(1,2) = tensor.yz();
  tensor_CMM(2,0) = tensor.xz();
  tensor_CMM(2,1) = tensor.yz();
  tensor_CMM(2,2) = tensor.zz();

  ColumnMajorMatrix tmp = _rot_mat * tensor_CMM;
  ColumnMajorMatrix rotT = _rot_mat.transpose();
  ColumnMajorMatrix rotated_tensor = tmp * rotT;

  return rotated_tensor;
}

ColumnMajorMatrix
CrackFrontRTheta::rotateToCrackFrontCoords(const ColumnMajorMatrix tensor) const
{
  ColumnMajorMatrix tmp = _rot_mat * tensor;
  ColumnMajorMatrix rotT = _rot_mat.transpose();
  ColumnMajorMatrix rotated_tensor = tmp * rotT;

  return rotated_tensor;
}

void 
CrackFrontRTheta::calcRTheta(Point & p, Point & crack_front_point, Point & crack_direction)
{
  RealVectorValue tangent_direction;
  tangent_direction(2) = 1.0;
  _crack_plane_normal = tangent_direction.cross(crack_direction);
  _rot_mat(0,0) = crack_direction(0);
  _rot_mat(0,1) = crack_direction(1);
  _rot_mat(0,2) = crack_direction(2);
  _rot_mat(1,0) = _crack_plane_normal(0);
  _rot_mat(1,1) = _crack_plane_normal(1);
  _rot_mat(1,2) = _crack_plane_normal(2);
  _rot_mat(2,0) = 0.0;
  _rot_mat(2,1) = 0.0;
  _rot_mat(2,2) = 0.0;
  _rot_mat(2,2) = 1.0;

  Point closest_point(0.0);
  RealVectorValue crack_front_point_rot = rotateToCrackFrontCoords(crack_front_point);

  RealVectorValue crack_front_edge = rotateToCrackFrontCoords(tangent_direction);

  Point p_rot = rotateToCrackFrontCoords(p);
  p_rot = p_rot - crack_front_point_rot;

  RealVectorValue closest_point_to_p = p_rot;

  //Find r, the distance between the qp and the crack front
  RealVectorValue r_vec = p_rot;
  _r = r_vec.size();

  //Find theta, the angle between r and the crack front plane
  RealVectorValue crack_plane_normal = rotateToCrackFrontCoords(_crack_plane_normal);
  Real p_to_plane_dist = std::abs(closest_point_to_p*crack_plane_normal);

  //Determine if p is above or below the crack plane
  Real y_local = p_rot(1) - closest_point(1);

  //Determine if p is in front of or behind the crack front
  RealVectorValue p2(p_rot);
  p2(1) = 0;
  RealVectorValue p2_vec = p2 - closest_point;
  Real ahead = crack_front_edge(2) * p2_vec(0) - crack_front_edge(0) * p2_vec(2);

  Real x_local(0);
  if (ahead >= 0)
    x_local = 1;
  else
    x_local = -1;

  //Calculate theta based on in which quadrant in the crack front coordinate
  //system the qp is located
  if (_r > 0)
  {
    Real theta_quadrant1(0.0);
    if (MooseUtils::absoluteFuzzyEqual(_r, p_to_plane_dist, 1e-10))
      theta_quadrant1 = 0.5*libMesh::pi;
    else if (p_to_plane_dist > _r)
      mooseError("Invalid distance p_to_plane_dist in CrackFrontDefinition::calculateRThetaToCrackFront");
    else
      theta_quadrant1 = std::asin(p_to_plane_dist/_r);

    if (x_local >= 0 && y_local >= 0)
      _theta = theta_quadrant1;

    else if (x_local < 0 && y_local >= 0)
      _theta = libMesh::pi - theta_quadrant1;

    else if (x_local < 0 && y_local < 0)
      _theta = -(libMesh::pi - theta_quadrant1);

    else if (x_local >= 0 && y_local < 0)
      _theta = -theta_quadrant1;
  }
  else if (_r == 0)
    _theta = 0;
  else
    mooseError("Invalid distance r in CrackFrontRTheta::calculateRTheta");

}

