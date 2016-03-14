/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "XFEMCircleCut2D.h"

#include "MooseError.h"
#include "libmesh/string_to_enum.h"

XFEMCircleCut2D::XFEMCircleCut2D(Real xc, Real yc, Real rc) :
    XFEMGeometricCut(0.0, 0.0),
    _xc(xc),
    _yc(yc),
    _rc(rc)
{
}

XFEMCircleCut2D::~XFEMCircleCut2D()
{
}

bool
XFEMCircleCut2D::active(Real time)
{
  return true;
}


bool
XFEMCircleCut2D::cutElementByGeometry(const Elem* elem, std::vector<CutEdge> & cut_edges, Real time)
{
  bool cut_elem = false;

  unsigned int n_sides = elem->n_sides();

  for (unsigned int i = 0; i < n_sides; ++i)
  {
    // This returns the lowest-order type of side, which should always
    // be an EDGE2 here because this class is for 2D only.
    UniquePtr<Elem> curr_side = elem->side(i);
    if (curr_side->type() != EDGE2)
      mooseError("In cutElementByGeometry element side must be EDGE2, but type is: " << libMesh::Utility::enum_to_string(curr_side->type())
                 << " base element type is: " << libMesh::Utility::enum_to_string(elem->type()));
    const Node *node1 = curr_side->get_node(0);
    const Node *node2 = curr_side->get_node(1);

    Real seg_originx = (*node1)(0);
    Real seg_originy = (*node1)(1);
    Real seg_endx = (*node2)(0);
    Real seg_endy = (*node2)(1);

    Real ls_origin = std::sqrt( (seg_originx - _xc)*(seg_originx - _xc) + (seg_originy - _yc)*(seg_originy - _yc) ) - _rc;	        
    Real ls_end = std::sqrt( (seg_endx - _xc)*(seg_endx - _xc) + (seg_endy - _yc)*(seg_endy - _yc) ) - _rc;

    if (ls_origin*ls_end<=-1.0e-10) //TODO: revisit end cases for intersections with corners
    {
      cut_elem = true;
      CutEdge mycut;
      mycut.id1 = node1->id();
      mycut.id2 = node2->id();
      Real seg_int_frac = std::abs(ls_origin)/std::abs(ls_origin - ls_end);
      mycut.distance = seg_int_frac;
      mycut.host_side_id = i;
      cut_edges.push_back(mycut);
    }
  }

  return cut_elem;
}

bool
XFEMCircleCut2D::cutElementByGeometry(const Elem* /*elem*/, std::vector<CutFace> & /*cut_faces*/, Real /*time*/)
{
  mooseError("invalid method for 2D mesh cutting");
  return false;
}

bool
XFEMCircleCut2D::cutFragmentByGeometry(std::vector<std::vector<Point> > & frag_edges,
                                       std::vector<CutEdge> & cut_edges, Real time)
{
  //Does not allow for double-cut 
  return false;
}

bool
XFEMCircleCut2D::cutFragmentByGeometry(std::vector<std::vector<Point> > & /*frag_faces*/,
                                       std::vector<CutFace> & /*cut_faces*/, Real /*time*/)
{
  mooseError("invalid method for 2D mesh cutting");
  return false;
}
