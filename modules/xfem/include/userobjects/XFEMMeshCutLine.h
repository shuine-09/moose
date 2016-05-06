/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef XFEMMESHCUTLINE_H
#define XFEMMESHCUTLINE_H

#include "XFEMMeshCutUserObject.h"

class XFEMMeshCutLine : public XFEMMeshCutUserObject
{
public:

  XFEMMeshCutLine(const InputParameters & parameters);

  virtual ~XFEMMeshCutLine() {}

protected:
  virtual bool cutElement(std::vector<CutEdge> & cut_edges);
  virtual bool cutElement(std::vector<CutFace> & cut_faces);
  virtual bool cutFragment(std::vector<std::vector<Point> > & frag_edges,
                           std::vector<CutEdge> & cut_edges);
  virtual bool cutFragment(std::vector<std::vector<Point> > & frag_faces,
                           std::vector<CutFace> & cut_faces);


  bool IntersectSegmentWithCutLine(const Point & segment_point1,
                                   const Point & segment_point2,
                                   const Point & cutting_line_point1,
                                   const Point & cutting_line_point2,
                                   Real & segment_intersection_fraction);
  Real crossProduct2D(Real ax, Real ay, Real bx, Real by);
private:
  std::vector<Real> _cut_data;
  Point _cut_line_start;
  Point _cut_line_end;
};

#endif //XFEMMESHCUTLINE_H
