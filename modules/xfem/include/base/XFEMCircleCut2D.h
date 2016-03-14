/********************:********************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef XFEMCIRCLECUT2D_H
#define XFEMCIRCLECUT2D_H

#include "XFEMGeometricCut.h"

class XFEMCircleCut2D : public XFEMGeometricCut
{
public:

  XFEMCircleCut2D(Real xc, Real yc, Real rc);
  ~XFEMCircleCut2D();

  virtual bool active(Real time);

  virtual bool cutElementByGeometry(const Elem* elem, std::vector<CutEdge> & cut_edges, Real time);
  virtual bool cutElementByGeometry(const Elem* elem, std::vector<CutFace> & cut_faces, Real time);

  virtual bool cutFragmentByGeometry(std::vector<std::vector<Point> > & frag_edges,
                                     std::vector<CutEdge> & cut_edges, Real time);
  virtual bool cutFragmentByGeometry(std::vector<std::vector<Point> > & frag_faces,
                                     std::vector<CutFace> & cut_faces, Real time);

private:

  Real _xc;
  Real _yc;
  Real _rc;
};

#endif
