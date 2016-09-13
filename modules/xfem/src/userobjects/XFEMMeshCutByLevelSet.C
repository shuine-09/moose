/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "XFEMMeshCutByLevelSet.h"
#include "MooseError.h"
#include "libmesh/string_to_enum.h"

template<>
InputParameters validParams<XFEMMeshCutByLevelSet>()
{
  InputParameters params = validParams<XFEMMeshCutUserObject>();
  params.addRequiredParam<NonlinearVariableName>("level_set_var", "The name of level set variable used to represent the interface");
  params.addClassDescription("XFEM mesh cut by level set function");
  return params;
}

XFEMMeshCutByLevelSet::XFEMMeshCutByLevelSet(const InputParameters & parameters) :
    XFEMMeshCutUserObject(parameters),
    _ls_var_name(parameters.get<NonlinearVariableName>("level_set_var")),
    _ls_var(_fe_problem.getVariable(_tid, _ls_var_name)),
    _aux_system(_fe_problem.getAuxiliarySystem()),
    _aux_solution(_aux_system.currentSolution())
{
}

void
XFEMMeshCutByLevelSet::execute()
{
  if (_current_elem->dim() == 2)
  {
    unsigned int _current_eid = _current_elem->id();
    std::map<unsigned int, RealVectorValue>::iterator mit1, mit2;
    mit1 = _marked_elems_host_id.find(_current_eid);
    mit2 = _marked_elems_distance.find(_current_eid);

    std::vector<CutEdge> cut_edges;
    if (cutElementByGeometry(_current_elem, cut_edges))
    {
      if (mit1 != _marked_elems_host_id.end())
      {
        mooseError("ERROR: element "<<_current_eid<<" already marked for cut.");
      }

      RealVectorValue distance(0.0, 0.0, 0.0);
      RealVectorValue host_id(0.0, 0.0, 0.0);

      for (unsigned int i = 0; i < cut_edges.size(); ++i)
      {
        distance(i) = cut_edges[i].distance;
        host_id(i) = cut_edges[i].host_side_id;
      }
      _marked_elems_distance[_current_eid] = distance;
      _marked_elems_host_id[_current_eid] = host_id;
    }
  }
  else
  {
    unsigned int _current_eid = _current_elem->id();

    std::vector<CutFace> cut_faces;
    if (cutElementByGeometry(_current_elem, cut_faces))
    {
      for (unsigned int i = 0; i < cut_faces.size(); ++i)
      {
        unsigned int k = cut_faces[i].face_id;

        RealVectorValue distance(0.0, 0.0, 0.0);
        RealVectorValue host_id(0.0, 0.0, 0.0);

        for (unsigned int l = 0; l < cut_faces[i].face_edge.size(); ++l)
        {
          host_id(l) = cut_faces[i].face_edge[l];
          distance(l) = cut_faces[i].position[l];
        }
        switch(k) {
          case 0 :
            {
              _marked_elems_distance[_current_eid] = distance;
              _marked_elems_host_id[_current_eid] = host_id;
              break;
            }
          case 1 :
            {
              _marked_elems_distance2[_current_eid] = distance;
              _marked_elems_host_id2[_current_eid] = host_id;
              break;
            }
          case 2 :
            {
              _marked_elems_distance3[_current_eid] = distance;
              _marked_elems_host_id3[_current_eid] = host_id;
              break;
            }
          case 3 :
            {
              _marked_elems_distance4[_current_eid] = distance;
              _marked_elems_host_id4[_current_eid] = host_id;
              break;
            }
          case 4 :
            {
              _marked_elems_distance5[_current_eid] = distance;
              _marked_elems_host_id5[_current_eid] = host_id;
              break;
            }
          case 5 :
            {
              _marked_elems_distance6[_current_eid] = distance;
              _marked_elems_host_id6[_current_eid] = host_id;
              break;
            }

          default : mooseError("ERROR: bad switch value in XFEMMEshCutByLevelSet");
        }
      }
    }
  }
}

bool
XFEMMeshCutByLevelSet::cutElementByGeometry(const Elem* elem, std::vector<CutFace> & cut_faces)
{
  bool cut_elem = false;

  for (unsigned int i = 0; i < elem->n_sides(); ++i)
  {
    // This returns the lowest-order type of side.
    std::unique_ptr<Elem> curr_side = elem->side(i);
    if (curr_side->dim() != 2)
      mooseError("In XFEMMeshCutByLevelSet dimension of side must be 2, but it is " << curr_side->dim());
    unsigned int n_edges = curr_side->n_sides();

    std::vector<unsigned int> cut_edges;
    std::vector<Real> cut_pos;

    for (unsigned int j = 0; j < n_edges; j++)
    {
      // This returns the lowest-order type of side.
      std::unique_ptr<Elem> curr_edge = curr_side->side(j);
      if (curr_edge->type() != EDGE2)
        mooseError("In XFEMMeshCutByLevelSet face edge must be EDGE2, but type is: " << libMesh::Utility::enum_to_string(curr_edge->type())
                   << " base element type is: " << libMesh::Utility::enum_to_string(elem->type()));
      
      const Node *node1 = curr_edge->get_node(0);
      const Node *node2 = curr_edge->get_node(1);

      dof_id_type ls_dof_id_1 = node1->dof_number(_aux_system.number(), _ls_var.number(), 0);
      dof_id_type ls_dof_id_2 = node2->dof_number(_aux_system.number(), _ls_var.number(), 0);

      Number ls_node_1 = (*_aux_solution)(ls_dof_id_1);
      Number ls_node_2 = (*_aux_solution)(ls_dof_id_2);

      if (ls_node_1 * ls_node_2 < 0)
      {
        Real seg_int_frac = std::abs(ls_node_1)/std::abs(ls_node_1 - ls_node_2);
        cut_edges.push_back(j);
        cut_pos.push_back(seg_int_frac);
      }
    }

    if (cut_edges.size() == 2)
    {
      cut_elem = true;
      CutFace mycut;
      mycut.face_id = i;
      mycut.face_edge.push_back(cut_edges[0]);
      mycut.face_edge.push_back(cut_edges[1]);
      mycut.position.push_back(cut_pos[0]);
      mycut.position.push_back(cut_pos[1]);
      cut_faces.push_back(mycut);
    }
  }
  return cut_elem;
}

bool
XFEMMeshCutByLevelSet::cutElementByGeometry(const Elem* elem, std::vector<CutEdge> & cut_edges)
{
  bool cut_elem = false;

  unsigned int n_sides = elem->n_sides();

  for (unsigned int i = 0; i < n_sides; ++i)
  {
    UniquePtr<Elem> curr_side = elem->side(i);

    if (curr_side->type() != EDGE2)
      mooseError("In XFEMMeshCutByLevelSet element side must be EDGE2, but type is: " << libMesh::Utility::enum_to_string(curr_side->type())
                 << " base element type is: " << libMesh::Utility::enum_to_string(elem->type()));

    const Node *node1 = curr_side->get_node(0);
    const Node *node2 = curr_side->get_node(1);

    dof_id_type ls_dof_id_1 = node1->dof_number(_aux_system.number(), _ls_var.number(), 0);
    dof_id_type ls_dof_id_2 = node2->dof_number(_aux_system.number(), _ls_var.number(), 0);

    Number ls_node_1 = (*_aux_solution)(ls_dof_id_1);
    Number ls_node_2 = (*_aux_solution)(ls_dof_id_2);

    if (ls_node_1 * ls_node_2 < 0)
    {
      cut_elem = true;
      CutEdge mycut;
      mycut.id1 = node1->id();
      mycut.id2 = node2->id();
      Real seg_int_frac = std::abs(ls_node_1)/std::abs(ls_node_1 - ls_node_2);
      mycut.distance = seg_int_frac;
      mycut.host_side_id = i;
      cut_edges.push_back(mycut);
    }
  }
  return cut_elem;
}

bool
XFEMMeshCutByLevelSet::cutFragmentByGeometry(std::vector<std::vector<Point> > & frag_edges,
                                          std::vector<CutEdge> & cut_edges)
{
  bool cut_frag = false;

  return cut_frag;
}

bool
XFEMMeshCutByLevelSet::cutFragmentByGeometry(std::vector<std::vector<Point> > & /*frag_faces*/,
                                          std::vector<CutFace> & /*cut_faces*/)
{
  mooseError("invalid method for 2D mesh cutting");
  return false;
}

