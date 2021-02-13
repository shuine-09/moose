//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"

// Forward Declarations
class InterfaceMeshCut3DUserObject;

/**
 * Calculate level set values for an interface that is defined by a set of line segments
 */
class MeshCutLevelSetAux : public AuxKernel
{
public:
  static InputParameters validParams();

  MeshCutLevelSetAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  /**
   * calculate the signed distance value for a given point.
   * @param p Coordinate of point
   * @return Signed distance
   */
  Real calculateSignedDistance(Point p);

  Real pointSegmentDistance(const Point & x0, const Point & x1, const Point & x2, Point & xp);
  Real pointTriangleDistance(const Point & x0,
                             const Point & x1,
                             const Point & x2,
                             const Point & x3,
                             Point & xp,
                             unsigned int & location_index);
  /// Pointer to the InterfaceMeshCut3DUserObject object
  const InterfaceMeshCut3DUserObject * _mesh_cut_uo;

  /// The structural mesh
  MooseMesh & _mesh;

  const std::map<unsigned int, std::array<Point, 7>> * _pseudo_normal;
};
