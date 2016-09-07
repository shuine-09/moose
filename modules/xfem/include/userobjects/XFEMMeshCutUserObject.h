/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef XFEMMESHCUTUSEROBJECT_H
#define XFEMMESHCUTUSEROBJECT_H

#include "DiscreteElementUserObject.h"
#include "XFEMGeometricCut.h"

class XFEMMeshCutUserObject;

template<>
InputParameters validParams<XFEMMeshCutUserObject>();

/**
 * XFEM mesh cut userobject base class.
 */
class XFEMMeshCutUserObject : public DiscreteElementUserObject
{
public:
  XFEMMeshCutUserObject(const InputParameters & parameters);

  virtual bool cutElementByGeometry(const Elem* elem, std::vector<CutEdge> & cut_edges) = 0;
  virtual bool cutElementByGeometry(const Elem* elem, std::vector<CutFace> & cut_faces) = 0;

  virtual bool cutFragmentByGeometry(std::vector<std::vector<Point> > & frag_edges,
                                     std::vector<CutEdge> & cut_edges) = 0;
  virtual bool cutFragmentByGeometry(std::vector<std::vector<Point> > & frag_faces,
                                     std::vector<CutFace> & cut_faces) = 0;
};

#endif // XFEMMESHCUTUSEROBJECT_H
