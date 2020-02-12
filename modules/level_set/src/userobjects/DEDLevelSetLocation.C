//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DEDLevelSetLocation.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "Function.h"
#include "Assembly.h"

#include "libmesh/quadrature.h"

registerMooseObject("TensorMechanicsApp", DEDLevelSetLocation);

defineLegacyParams(DEDLevelSetLocation);

InputParameters
DEDLevelSetLocation::validParams()
{
  InputParameters params = ElementUserObject::validParams();
  params.addClassDescription("Find laser spot.");
  params.addRequiredCoupledVar("level_set", "Level set variable");
  params.addRequiredParam<VariableName>(
      "level_set_var", "The name of level set variable used to represent the interface");
  params.addParam<FunctionName>("laser_center_x", 0, "The laser center function of x coordinate.");
  params.addParam<FunctionName>("laser_center_y", 0, "The laser center function of y coordinate.");
  params.addParam<FunctionName>("laser_center_z", 0, "The laser center function of z coordinate.");

  return params;
}

DEDLevelSetLocation::DEDLevelSetLocation(const InputParameters & parameters)
  : ElementUserObject(parameters),
    _laser_center_x(getFunction("laser_center_x")),
    _laser_center_y(getFunction("laser_center_y")),
    _laser_center_z(getFunction("laser_center_z")),
    _level_set_var_number(_subproblem
                              .getVariable(_tid,
                                           parameters.get<VariableName>("level_set_var"),
                                           Moose::VarKindType::VAR_ANY,
                                           Moose::VarFieldType::VAR_FIELD_STANDARD)
                              .number()),
    _system(_subproblem.getSystem(getParam<VariableName>("level_set_var"))),
    _solution(_system.current_local_solution.get())
{
}

void
DEDLevelSetLocation::initialize()
{
  _location_x = 0;
  _location_y = 0;
  _location_z = 0;
}

void
DEDLevelSetLocation::execute()
{
  bool cut_elem = false;
  bool at_laser = false;

  Real x = _laser_center_x.value(_t, _current_elem->centroid());

  unsigned int n_sides = _current_elem->n_sides();
  Real seg_int_frac = 0;

  for (unsigned int i = 0; i < n_sides; ++i)
  {
    UniquePtr<const Elem> curr_side = _current_elem->side_ptr(i);

    const Node * node1 = curr_side->node_ptr(0);
    const Node * node2 = curr_side->node_ptr(1);

    dof_id_type ls_dof_id_1 = node1->dof_number(_system.number(), _level_set_var_number, 0);
    dof_id_type ls_dof_id_2 = node2->dof_number(_system.number(), _level_set_var_number, 0);

    Number ls_node_1 = (*_solution)(ls_dof_id_1)-0.5;
    Number ls_node_2 = (*_solution)(ls_dof_id_2)-0.5;

    if (ls_node_1 * ls_node_2 <= 0)
    {
      cut_elem = true;
      // std::cout << "node1 = " << *node1 << ", node2 = " << *node2 << std::endl;
      seg_int_frac = std::abs(ls_node_1) / std::abs(ls_node_1 - ls_node_2 + 1.0e-10) *
                         std::abs((*node1)(1) - (*node2)(1)) +
                     (*node1)(1);
    }
  }

  for (unsigned int i = 0; i < n_sides; ++i)
  {
    UniquePtr<const Elem> curr_side = _current_elem->side_ptr(i);

    const Node * node1 = curr_side->node_ptr(0);
    const Node * node2 = curr_side->node_ptr(1);

    Real left = std::min((*node1)(0), (*node2)(0));
    Real right = std::max((*node1)(0), (*node2)(0));

    if (x >= left && x < right)
    {
      at_laser = true;
      // std::cout << "at_laser left " << left << ", right = " << right << std::endl;
    }
  }

  if (at_laser && cut_elem && _location_y == 0 && _location_x == 0)
  {
    _location_x += x;
    _location_y += seg_int_frac; //(_current_elem->centroid())(1);
    _location_z += 0;
  }
  else
  {
    _location_x += 0;
    _location_y += 0;
    _location_z += 0;
  }
}

void
DEDLevelSetLocation::threadJoin(const UserObject & uo)
{
  const DEDLevelSetLocation & ls_uo = static_cast<const DEDLevelSetLocation &>(uo);

  _location_x += ls_uo._location_x;
  _location_y += ls_uo._location_y;
  _location_z += ls_uo._location_z;
}

void
DEDLevelSetLocation::finalize()
{
  gatherSum(_location_x);
  gatherSum(_location_y);
  gatherSum(_location_z);

  std::cout << "Location = " << getLaserSpotLocation() << std::endl;
}
