/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef XFEMMESHCUTUSEROBJECT_H
#define XFEMMESHCUTUSEROBJECT_H

#include "ElementUserObject.h"
#include "XFEMGeometricCut.h"

class XFEM;

class XFEMMeshCutUserObject;

template<>
InputParameters validParams<XFEMMeshCutUserObject>();

/**
 * XFEM mesh cut userobject base class.
 */
class XFEMMeshCutUserObject : public ElementUserObject
{
public:
  XFEMMeshCutUserObject(const InputParameters & parameters);

  virtual void initialize();
  virtual void execute();
  virtual void threadJoin(const UserObject &y);
  virtual void finalize();


  virtual bool cutElementByGeometry(const Elem* elem, std::vector<CutEdge> & cut_edges) = 0;
  virtual bool cutElementByGeometry(const Elem* elem, std::vector<CutFace> & cut_faces) = 0;

  virtual bool cutFragmentByGeometry(std::vector<std::vector<Point> > & frag_edges,
                                     std::vector<CutEdge> & cut_edges) = 0;
  virtual bool cutFragmentByGeometry(std::vector<std::vector<Point> > & frag_faces,
                                     std::vector<CutFace> & cut_faces) = 0;

private:
  MooseMesh & _mesh;
  MooseSharedPointer<XFEM> _xfem;
protected:
  std::map<unsigned int, RealVectorValue > _marked_elems_distance;
  std::map<unsigned int, RealVectorValue > _marked_elems_host_id;

};

#endif // XFEMMESHCUTUSEROBJECT_H
