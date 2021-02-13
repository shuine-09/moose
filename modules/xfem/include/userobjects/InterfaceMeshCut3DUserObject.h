//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeometricCutUserObject.h"
#include "XFEMMovingInterfaceVelocityBase.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/exodusII_io_helper.h"

#include <array>

class Function;
// Forward declarations
class PointValueAtXFEMInterface;

/**
 * InterfaceMeshCut3DUserObject: (1) reads in a mesh describing the crack surface,
 * (2) uses the mesh to do initial cutting of 3D elements, and
 * (3) grows the mesh based on prescribed growth functions.
 */

class InterfaceMeshCut3DUserObject : public GeometricCutUserObject
{
public:
  static InputParameters validParams();

  InterfaceMeshCut3DUserObject(const InputParameters & parameters);

  virtual void initialSetup() override;
  virtual void initialize() override;
  virtual const std::vector<Point>
  getCrackFrontPoints(unsigned int num_crack_front_points) const override;

  std::shared_ptr<MeshBase> getCutMesh() const { return _cut_mesh; };
  const auto * getPseudoNormal() const { return &_pseudo_normal; };

  virtual bool cutElementByGeometry(const Elem * elem,
                                    std::vector<Xfem::CutEdge> & cut_edges,
                                    std::vector<Xfem::CutNode> & cut_nodes,
                                    Real time) const override;
  virtual bool cutElementByGeometry(const Elem * elem,
                                    std::vector<Xfem::CutFace> & cut_faces,
                                    Real time) const override;
  virtual bool cutFragmentByGeometry(std::vector<std::vector<Point>> & frag_edges,
                                     std::vector<Xfem::CutEdge> & cut_edges,
                                     Real time) const override;
  virtual bool cutFragmentByGeometry(std::vector<std::vector<Point>> & frag_faces,
                                     std::vector<Xfem::CutFace> & cut_faces,
                                     Real time) const override;

protected:
  /// The cutter mesh
  // std::unique_ptr<MeshBase> _cut_mesh;
  std::shared_ptr<MeshBase> _cut_mesh;

  std::shared_ptr<ExodusII_IO> _exodus_io;

  /// The cutter mesh must be 2D
  const unsigned int _cut_elem_dim = 2;

  /// The structural mesh
  MooseMesh & _mesh;

  /// The structural mesh must be 3D only
  const unsigned int _elem_dim = 3;

  /**
    Check if a line intersects with an element
   */
  virtual bool intersectWithEdge(const Point & p1,
                                 const Point & p2,
                                 const std::vector<Point> & _vertices,
                                 Point & pint) const;

  Point findNormalatNode(const Node & node);

  /**
    Check if point p is inside the edge p1-p2
   */
  bool isInsideEdge(const Point & p1, const Point & p2, const Point & p) const;

  /**
    Get the relative position of p from p1
   */
  Real getRelativePosition(const Point & p1, const Point & p2, const Point & p) const;

  /**
    Check if point p is inside a plane
   */
  bool isInsideCutPlane(const std::vector<Point> & _vertices, const Point & p) const;

  std::map<unsigned int, std::array<Point, 7>> _pseudo_normal;

  /// Pointer to XFEMMovingInterfaceVelocityBase object
  const XFEMMovingInterfaceVelocityBase * _interface_velocity;
};
