/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "LevelSetCutUserObject.h"
#include "libmesh/string_to_enum.h"
#include "SystemBase.h"
#include "MooseError.h"
#include "MooseVariable.h"

template <>
InputParameters
validParams<LevelSetCutUserObject>()
{
  // Get input parameters from parent class
  InputParameters params = validParams<GeometricCut2DUserObject>();

  // Add required parameters
  params.addRequiredParam<AuxVariableName>(
      "level_set_var", "The name of level set variable used to represent the interface");
  // Class description
  params.addClassDescription("Creates a UserObject for a line segment cut on 2D meshes for XFEM");
  // Return the parameters
  return params;
}

LevelSetCutUserObject::LevelSetCutUserObject(const InputParameters & parameters)
  : GeometricCut2DUserObject(parameters),
    _ls_var_name(parameters.get<AuxVariableName>("level_set_var")),
    _aux_system(_fe_problem.getAuxiliarySystem()),
    _ls_var(_fe_problem.getVariable(_tid, _ls_var_name)),
    _aux_solution(_aux_system.currentSolution())
{
}

bool
LevelSetCutUserObject::cutElementByGeometry(const Elem * elem,
                                            std::vector<CutEdge> & cut_edges,
                                            Real time) const
{
  NumericVector<Number> & current_aux_solution = _aux_system.serializedSolution();
  _aux_system.serializeSolution();

  bool cut_elem = false;

  unsigned int n_sides = elem->n_sides();

  for (unsigned int i = 0; i < n_sides; ++i)
  {
    UniquePtr<Elem> curr_side = elem->side(i);

    // if (curr_side->type() != EDGE2)
    //   mooseError("In XFEMMeshCutByLevelSet element side must be EDGE2, but type is: "
    //              << libMesh::Utility::enum_to_string(curr_side->type())
    //              << " base element type is: "
    //              << libMesh::Utility::enum_to_string(elem->type()));

    const Node * node1 = curr_side->get_node(0);
    const Node * node2 = curr_side->get_node(1);

    dof_id_type ls_dof_id_1 = node1->dof_number(_aux_system.number(), _ls_var.number(), 0);
    dof_id_type ls_dof_id_2 = node2->dof_number(_aux_system.number(), _ls_var.number(), 0);

    Number ls_node_1 = (current_aux_solution)(ls_dof_id_1);
    Number ls_node_2 = (current_aux_solution)(ls_dof_id_2);

    if (ls_node_1 * ls_node_2 < 0)
    {
      cut_elem = true;
      CutEdge mycut;
      mycut.id1 = node1->id();
      mycut.id2 = node2->id();
      Real seg_int_frac = std::abs(ls_node_1) / std::abs(ls_node_1 - ls_node_2);
      mycut.distance = seg_int_frac;
      mycut.host_side_id = i;
      cut_edges.push_back(mycut);
    }
  }
  current_aux_solution.close();
  return cut_elem;
}
