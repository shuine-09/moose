//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MeshCutLevelSetAux.h"
#include "InterfaceMeshCut3DUserObject.h"
#include "XFEMFuncs.h"
#include "MooseError.h"
#include "libmesh/string_to_enum.h"
#include "MooseMesh.h"
#include "libmesh/face_tri3.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/plane.h"
#include "Function.h"

registerMooseObject("XFEMApp", MeshCutLevelSetAux);

InputParameters
MeshCutLevelSetAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription(
      "Auxiliary Kernel that calcuates level set value using line segments' description.");
  params.addParam<UserObjectName>(
      "mesh_cut_user_object", "Name of GeometricCutUserObject that gives cut mesh information.");
  return params;
}

MeshCutLevelSetAux::MeshCutLevelSetAux(const InputParameters & parameters)
  : AuxKernel(parameters), _mesh(_subproblem.mesh())
{
  if (!isNodal())
    mooseError("MeshCutLevelSetAux: Aux variable must be nodal variable.");

  FEProblemBase * fe_problem = dynamic_cast<FEProblemBase *>(&_subproblem);

  const UserObject * uo =
      &(fe_problem->getUserObjectBase(getParam<UserObjectName>("mesh_cut_user_object")));

  if (dynamic_cast<const InterfaceMeshCut3DUserObject *>(uo) == nullptr)
    mooseError("Failed to cast UserObject to InterfaceMeshCut3DUserObject in MeshCutLevelSetAux");

  _mesh_cut_uo = dynamic_cast<const InterfaceMeshCut3DUserObject *>(uo);
  _pseudo_normal = _mesh_cut_uo->getPseudoNormal();
}

Real
MeshCutLevelSetAux::pointSegmentDistance(const Point & x0,
                                         const Point & x1,
                                         const Point & x2,
                                         Point & xp)
{
  Point dx = x2 - x1;
  Real m2 = dx * dx;
  // find parameter value of closest point on segment
  Real s12 = (x2 - x0) * dx / m2;
  if (s12 < 0)
  {
    s12 = 0;
  }
  else if (s12 > 1)
  {
    s12 = 1;
  }
  // and find the distance
  xp = s12 * x1 + (1 - s12) * x2;
  return std::sqrt((x0 - xp) * (x0 - xp));
}

Real
MeshCutLevelSetAux::pointTriangleDistance(const Point & x0,
                                          const Point & x1,
                                          const Point & x2,
                                          const Point & x3,
                                          Point & xp,
                                          unsigned int & location_index)
{ // first find barycentric coordinates of closest point on infinite plane
  Point x13 = x1 - x3, x23 = x2 - x3, x03 = x0 - x3;
  Real m13 = x13 * x13, m23 = x23 * x23, d = x13 * x23;
  Real invdet = 1.0 / std::max(m13 * m23 - d * d, 1e-30);
  Real a = x13 * x03, b = x23 * x03;
  // the barycentric coordinates themselves
  Real w23 = invdet * (m23 * a - d * b);
  Real w31 = invdet * (m13 * b - d * a);
  Real w12 = 1 - w23 - w31;
  if (w23 >= 0 && w31 >= 0 && w12 >= 0)
  { // if we're inside the triangle
    location_index = 0;
    xp = w23 * x1 + w31 * x2 + w12 * x3;
    return std::sqrt((x0 - xp) * (x0 - xp));
  }
  else
  {              // we have to clamp to one of the edges
    if (w23 > 0) // this rules out edge 2-3 for us
    {
      Point xp1, xp2;
      Real distance_12 = pointSegmentDistance(x0, x1, x2, xp1);
      Real distance_13 = pointSegmentDistance(x0, x1, x3, xp2);
      Real distance_1 = std::sqrt((x0 - x1) * (x0 - x1));
      if (std::min(distance_12, distance_13) < distance_1)
      {
        if (distance_12 < distance_13)
        {
          location_index = 4;
          xp = xp1;
          return distance_12;
        }
        else
        {
          location_index = 6;
          xp = xp2;
          return distance_13;
        }
      }
      else
      {
        location_index = 1;
        xp = x1;
        return distance_1;
      }
    }
    else if (w31 > 0) // this rules out edge 1-3
    {
      Point xp1, xp2;
      Real distance_12 = pointSegmentDistance(x0, x1, x2, xp1);
      Real distance_23 = pointSegmentDistance(x0, x2, x3, xp2);
      Real distance_2 = std::sqrt((x0 - x2) * (x0 - x2));
      if (std::min(distance_12, distance_23) < distance_2)
      {
        if (distance_12 < distance_23)
        {
          location_index = 4;
          xp = xp1;
          return distance_12;
        }
        else
        {
          location_index = 5;
          xp = xp2;
          return distance_23;
        }
      }
      else
      {
        location_index = 2;
        xp = x2;
        return distance_2;
      }
    }
    else // w12 must be >0, ruling out edge 1-2
    {
      Point xp1, xp2;
      Real distance_23 = pointSegmentDistance(x0, x2, x3, xp1);
      Real distance_31 = pointSegmentDistance(x0, x3, x1, xp2);
      Real distance_3 = std::sqrt((x0 - x3) * (x0 - x3));
      if (std::min(distance_23, distance_31) < distance_3)
      {
        if (distance_23 < distance_31)
        {
          location_index = 5;
          xp = xp1;
          return distance_23;
        }
        else
        {
          location_index = 6;
          xp = xp2;
          return distance_31;
        }
      }
      else
      {
        location_index = 3;
        xp = x3;
        return distance_3;
      }
    }
  }
}

Real
MeshCutLevelSetAux::calculateSignedDistance(Point p)
{
  std::shared_ptr<MeshBase> cut_mesh = _mesh_cut_uo->getCutMesh();
  Real min_dist = std::numeric_limits<Real>::max();
  Real sgn = 1.0;
  for (const auto & cut_elem : cut_mesh->element_ptr_range())
  {
    std::vector<Point> vertices{
        cut_elem->node_ref(0), cut_elem->node_ref(1), cut_elem->node_ref(2)};
    unsigned int location_index;
    Point xp;
    Real dist = pointTriangleDistance(
        p, cut_elem->node_ref(0), cut_elem->node_ref(1), cut_elem->node_ref(2), xp, location_index);

    if (dist < std::abs(min_dist))
    {
      min_dist = dist;
      Point normal = ((*_pseudo_normal).find(cut_elem->id())->second)[location_index];
      if (normal * (p - xp) < 0.0)
        sgn = -1.0;
      else
        sgn = 1.0;
    }
  }
  return min_dist * sgn;
}

Real
MeshCutLevelSetAux::computeValue()
{
  return calculateSignedDistance(*_current_node);
}
