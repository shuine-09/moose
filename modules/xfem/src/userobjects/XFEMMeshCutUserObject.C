/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "XFEMMeshCutUserObject.h"

#include "XFEM.h"
#include "ElementFragmentAlgorithm.h"
#include "EFAElement.h"
#include "EFAElement2D.h"
#include "EFAElement3D.h"
#include "EFAFragment2D.h"
#include "EFAFragment3D.h"
#include "MooseMesh.h"

#include "libmesh/parallel_algebra.h"
#include "libmesh/parallel.h"

template<>
InputParameters validParams<XFEMMeshCutUserObject>()
{
  InputParameters params = validParams<ElementUserObject>();
  return params;
}

XFEMMeshCutUserObject::XFEMMeshCutUserObject(const InputParameters & parameters) :
    ElementUserObject(parameters),
    _mesh(_subproblem.mesh())
{
  FEProblem * fe_problem = dynamic_cast<FEProblem *>(&_subproblem);
  if (fe_problem == NULL)
    mooseError("Problem casting _subproblem to FEProblem in XFEMMeshCutUserObject");
  _xfem = MooseSharedNamespace::dynamic_pointer_cast<XFEM>(fe_problem->getXFEM());
  if (_xfem == NULL)
    mooseError("Problem casting to XFEM in XFEMMeshCutUserObject");
  if (isNodal())
    mooseError("XFEMMeshCutUserObject can only be run on an element variable");

  _efa_mesh = _xfem->getEFAMesh();
}

void
XFEMMeshCutUserObject::initialize()
{
}

bool
XFEMMeshCutUserObject::cutXFEMMesh2D()
{
  bool marked_edges = false;
  
  std::vector<CutEdge> elem_cut_edges;
  std::vector<CutEdge> frag_cut_edges;
  std::vector<std::vector<Point> > frag_edges;
  EFAElement * EFAelem = _efa_mesh->getElemByID(_current_elem->id());
  EFAElement2D * CEMElem = dynamic_cast<EFAElement2D*>(EFAelem);

  if (!CEMElem)
    mooseError("EFAelem is not of EFAelement2D type");

  // continue if elem has been already cut twice - IMPORTANT
  if (CEMElem->isFinalCut())
    return marked_edges;

  // get fragment edges
  _xfem->getFragmentEdges(_current_elem, CEMElem, frag_edges);

  cutElement(elem_cut_edges);

  if (CEMElem->numFragments() > 0)
    cutFragment(frag_edges, frag_cut_edges);

  for (unsigned int i = 0; i < elem_cut_edges.size(); ++i) // mark element edges
  {
    if (!CEMElem->isEdgePhantom(elem_cut_edges[i].host_side_id)) // must not be phantom edge
    {
      _efa_mesh->addElemEdgeIntersection(_current_elem->id(), elem_cut_edges[i].host_side_id,
                                        elem_cut_edges[i].distance);
      marked_edges = true;
    }
  }

  for (unsigned int i = 0; i < frag_cut_edges.size(); ++i) // MUST DO THIS AFTER MARKING ELEMENT EDGES
  {
    if (!CEMElem->getFragment(0)->isSecondaryInteriorEdge(frag_cut_edges[i].host_side_id))
    {
      _efa_mesh->addFragEdgeIntersection(_current_elem->id(), frag_cut_edges[i].host_side_id,
                                        frag_cut_edges[i].distance);
      marked_edges = true;
    }
  }
  return marked_edges;
}

bool
XFEMMeshCutUserObject::cutXFEMMesh3D()
{
  // TODO: 3D fragment cut
  bool marked_faces = false;
  
  std::vector<CutFace> elem_cut_faces;
  EFAElement * EFAelem = _efa_mesh->getElemByID(_current_elem->id());
  EFAElement3D * CEMElem = dynamic_cast<EFAElement3D*>(EFAelem);
  
  if (!CEMElem)
    mooseError("EFAelem is not of EFAelement3D type");

  // continue if elem has been already cut twice - IMPORTANT
  if (CEMElem->isFinalCut())
    return marked_faces;

  // mark cut faces for the element and its fragment
  cutElement(elem_cut_faces);

  for (unsigned int i = 0; i < elem_cut_faces.size(); ++i) // mark element faces
  {
    if (!CEMElem->isFacePhantom(elem_cut_faces[i].face_id)) // must not be phantom face
    {
      _efa_mesh->addElemFaceIntersection(_current_elem->id(), elem_cut_faces[i].face_id,
                                        elem_cut_faces[i].face_edge, elem_cut_faces[i].position);
      marked_faces = true;
    }
  }

  return marked_faces;
}

void
XFEMMeshCutUserObject::execute()
{
  if (_mesh.dimension() == 2)
    cutXFEMMesh2D();
  else if (_mesh.dimension() == 3)
    cutXFEMMesh3D();
}

void
XFEMMeshCutUserObject::threadJoin(const UserObject &y)
{
}

void
XFEMMeshCutUserObject::finalize()
{
}
