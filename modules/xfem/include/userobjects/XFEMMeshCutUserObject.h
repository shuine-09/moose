/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef XFEMMESHCUTUSEROBJECT_H
#define XFEMMESHCUTUSEROBJECT_H

#include "ElementUserObject.h"

class XFEM;
class ElementFragmentAlgorithm;

struct CutEdge
{
  unsigned int id1;
  unsigned int id2;
  Real distance;
  unsigned int host_side_id;
};

struct CutFace
{
  unsigned int face_id;
  std::vector<unsigned int> face_edge;
  std::vector<Real> position;
};

class XFEMMeshCutUserObject : public ElementUserObject
{
public:

  XFEMMeshCutUserObject(const InputParameters & parameters);

  virtual ~XFEMMeshCutUserObject() {}

  virtual void initialize();
  virtual void execute();
  virtual void threadJoin(const UserObject &y);
  virtual void finalize();

protected:
  virtual bool cutXFEMMesh2D();
  virtual bool cutXFEMMesh3D();
  virtual bool cutElement(std::vector<CutEdge> & cut_edges){return false;};
  virtual bool cutElement(std::vector<CutFace> & cut_faces){return false;};
  virtual bool cutFragment(std::vector<std::vector<Point> > & frag_edges,
                           std::vector<CutEdge> & cut_edges){return false;};
  virtual bool cutFragment(std::vector<std::vector<Point> > & frag_faces,
                           std::vector<CutFace> & cut_faces){return false;};

  MooseMesh & _mesh;
  MooseSharedPointer<XFEM> _xfem;
  ElementFragmentAlgorithm * _efa_mesh;
};

template<>
InputParameters validParams<XFEMMeshCutUserObject>();

#endif //XFEMMESHCUTUSEROBJECT_H
