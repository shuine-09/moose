/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "XFEMMeshCutLine.h"
#include "MooseError.h"
#include "libmesh/string_to_enum.h"

template<>
InputParameters validParams<XFEMMeshCutLine>()
{
  InputParameters params = validParams<XFEMMeshCutUserObject>();
  params.addParam<std::vector<Real> >("cut_data","Data for XFEM geometric cuts");
  return params;
}

XFEMMeshCutLine::XFEMMeshCutLine(const InputParameters & parameters) :
    XFEMMeshCutUserObject(parameters),
    _cut_data(getParam<std::vector<Real> >("cut_data"))
{
  Real x0 = _cut_data[0];
  Real y0 = _cut_data[1];
  Real x1 = _cut_data[2];
  Real y1 = _cut_data[3];
  _cut_line_start = Point(x0, y0, 0.0);
  _cut_line_end = Point(x1, y1, 0.0);
}

bool 
XFEMMeshCutLine::cutElement(std::vector<CutEdge> & cut_edges)
{
  bool cut_elem = false;

  unsigned int n_sides = _current_elem->n_sides();

  for (unsigned int i = 0; i < n_sides; ++i)
  {
    // This returns the lowest-order type of side, which should always
    // be an EDGE2 here because this class is for 2D only.
    UniquePtr<Elem> curr_side =_current_elem->side(i);
    if (curr_side->type() != EDGE2)
      mooseError("In cutElementByGeometry element side must be EDGE2, but type is: " << libMesh::Utility::enum_to_string(curr_side->type())
                 << " base element type is: " << libMesh::Utility::enum_to_string(_current_elem->type()));

    const Node *node1 = curr_side->get_node(0);
    const Node *node2 = curr_side->get_node(1);
    Real seg_int_frac = 0.0;

    if (IntersectSegmentWithCutLine(*node1, *node2, _cut_line_start, _cut_line_end, seg_int_frac))
    {
      cut_elem = true;
      CutEdge mycut;
      mycut.id1 = node1->id();
      mycut.id2 = node2->id();
      mycut.distance = seg_int_frac;
      mycut.host_side_id = i;
      cut_edges.push_back(mycut);
    }
  }
  return cut_elem;
}

bool 
XFEMMeshCutLine::cutFragment(std::vector<std::vector<Point> > & frag_edges, std::vector<CutEdge> & cut_edges)
{
  bool cut_frag = false;

  unsigned int n_sides = frag_edges.size();

  for (unsigned int i = 0; i < n_sides; ++i) 
  {                              
    Real seg_int_frac = 0.0;     

    if (IntersectSegmentWithCutLine(frag_edges[i][0], frag_edges[i][1], _cut_line_start, _cut_line_end, seg_int_frac))
    {                            
      cut_frag = true;
      CutEdge mycut;
      mycut.id1 = i;
      mycut.id2 = (i < (n_sides - 1) ? (i + 1) : 0);
      mycut.distance = seg_int_frac;
      mycut.host_side_id = i;
      cut_edges.push_back(mycut);
    }
  }

  return cut_frag;
}

bool
XFEMMeshCutLine::cutElement(std::vector<CutFace> & /*cut_faces*/)
{
  mooseError("invalid method for 2D mesh cutting");
  return false;
}

bool
XFEMMeshCutLine::cutFragment(std::vector<std::vector<Point> > & /*frag_faces*/, std::vector<CutFace> & /*cut_faces*/)
{
  mooseError("invalid method for 2D mesh cutting");
  return false;
}

bool
XFEMMeshCutLine::IntersectSegmentWithCutLine(const Point & segment_point1,
                                                const Point & segment_point2,
                                                const Point & cutting_line_point1,
                                                const Point & cutting_line_point2,
                                                Real & segment_intersection_fraction)
{
  // Use the algorithm described here to determine whether a line segment is intersected
  // by a cutting line, and to compute the fraction along that line where the intersection
  // occurs:
  // http://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect

  bool cut_segment = false;
  Real tol = 1.e-15;
  Point seg_dir = segment_point2 - segment_point1;
  Point cut_dir = cutting_line_point2 - cutting_line_point1;
  Point cut_start_to_seg_start = segment_point1 - cutting_line_point1;

  Real cut_dir_cross_seg_dir = crossProduct2D(cut_dir(0), cut_dir(1), seg_dir(0), seg_dir(1));

  if (std::abs(cut_dir_cross_seg_dir) > tol)
  {
    //Fraction of the distance along the cutting segment where it intersects the edge segment
    Real cut_int_frac = crossProduct2D(cut_start_to_seg_start(0), cut_start_to_seg_start(1), seg_dir(0), seg_dir(1)) /
      cut_dir_cross_seg_dir;

    if (cut_int_frac >= 0.0 && cut_int_frac <= 1.0)
    { //Cutting segment intersects the line of the edge segment, but the intersection point may be outside the segment
      Real int_frac = crossProduct2D(cut_start_to_seg_start(0), cut_start_to_seg_start(1), cut_dir(0), cut_dir(1)) /
        cut_dir_cross_seg_dir;
      if (int_frac >= 0.0 && int_frac <= 1.0) //TODO: revisit end cases for intersections with corners
      {
        cut_segment = true;
        segment_intersection_fraction = int_frac;
      }
    }
  }
  return cut_segment;
}

Real XFEMMeshCutLine::crossProduct2D(Real ax, Real ay, Real bx, Real by)
{
  return (ax * by - bx * ay);
} 
