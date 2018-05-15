//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PointValueAtXFEMInterface.h"

// MOOSE includes
#include "MooseMesh.h"
#include "MooseVariableFE.h"
#include "XFEM.h"

#include "libmesh/mesh_tools.h"
#include "MovingLineSegmentCutSetUserObject.h"

registerMooseObject("XFEMApp", PointValueAtXFEMInterface);

template <>
InputParameters
validParams<PointValueAtXFEMInterface>()
{
  InputParameters params = validParams<GeneralUserObject>();
  params.addRequiredCoupledVar(
      "variable", "The names of the variables that this VectorPostprocessor operates on");
  params.addParam<UserObjectName>(
      "geometric_cut_userobject",
      "Name of GeometricCutUserObject associated with this constraint.");
  params.addRequiredParam<VariableName>(
      "level_set_var", "The name of level set variable used to represent the interface");

  return params;
}

PointValueAtXFEMInterface::PointValueAtXFEMInterface(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    CoupleableMooseVariableDependencyIntermediateInterface(this, false),
    MooseVariableInterface<Real>(this, false),
    _mesh(_subproblem.mesh()),
    _level_set_var_number(
        _subproblem.getVariable(_tid, parameters.get<VariableName>("level_set_var")).number()),
    _system(_subproblem.getSystem(getParam<VariableName>("level_set_var"))),
    _solution(_system.current_local_solution.get())
{
  const UserObject * uo =
      &(_fe_problem.getUserObjectBase(getParam<UserObjectName>("geometric_cut_userobject")));

  if (dynamic_cast<const MovingLineSegmentCutSetUserObject *>(uo) == nullptr)
    mooseError("UserObject casting to GeometricCutUserObject in XFEMSingleVariableConstraint");

  _geo_cut = dynamic_cast<const MovingLineSegmentCutSetUserObject *>(uo);

  addMooseVariableDependency(mooseVariable());

  // TODO multiple variables
  std::vector<std::string> var_names(_coupled_moose_vars.size());

  for (unsigned int i = 0; i < _coupled_moose_vars.size(); i++)
    var_names[i] = _coupled_moose_vars[i]->name();

  _values_positive_level_set_side.resize(_points.size());
  _values_negative_level_set_side.resize(_points.size());
  _grad_values_positive_level_set_side.resize(_points.size());
  _grad_values_negative_level_set_side.resize(_points.size());

  _xfem = MooseSharedNamespace::dynamic_pointer_cast<XFEM>(_fe_problem.getXFEM());
  if (_xfem == nullptr)
    mooseError("Problem casting to XFEM in PointValueAtXFEMInterface");

  _elem_pairs = _xfem->getXFEMCutElemPairs(_xfem->getGeometricCutID(_geo_cut));

  std::vector<Real> cut_data = _geo_cut->getCutData();

  for (unsigned int i = 0; i < cut_data.size() / 6; ++i)
  {
    _points.push_back(Point(cut_data[i * 6 + 0], cut_data[i * 6 + 1]));
    if (i == cut_data.size() / 6 - 1)
      _points.push_back(Point(cut_data[i * 6 + 2], cut_data[i * 6 + 3]));
  }
}

void
PointValueAtXFEMInterface::initialize()
{
  // SamplerBase::initialize();

  // We do this here just in case it's been destroyed and recreated because of mesh adaptivity.
  _pl = _mesh.getPointLocator();

  // We may not find a requested point on a distributed mesh, and
  // that's okay.
  _pl->enable_out_of_mesh_mode();

  // Reset the point arrays
  _found_points.assign(_points.size(), false);
}

void
PointValueAtXFEMInterface::execute()
{
  _found_points.clear();
  _values_positive_level_set_side.clear();
  _values_negative_level_set_side.clear();
  _grad_values_positive_level_set_side.clear();
  _grad_values_negative_level_set_side.clear();
  _points.clear();

  std::vector<Real> cut_data = _geo_cut->getCutData();

  for (unsigned int i = 0; i < cut_data.size() / 6; ++i)
  {
    _points.push_back(Point(cut_data[i * 6 + 0], cut_data[i * 6 + 1]));
    if (i == cut_data.size() / 6 - 1)
      _points.push_back(Point(cut_data[i * 6 + 2], cut_data[i * 6 + 3]));
  }

  _found_points.assign(_points.size(), false);

  _values_positive_level_set_side.resize(_points.size());
  _values_negative_level_set_side.resize(_points.size());
  _grad_values_positive_level_set_side.resize(_points.size());
  _grad_values_negative_level_set_side.resize(_points.size());

  BoundingBox bbox = _mesh.getInflatedProcessorBoundingBox();

  /// So we don't have to create and destroy this
  std::vector<Point> point_vec(1);

  for (auto i = beginIndex(_points); i < _points.size(); ++i)
  {
    Point & p = _points[i];

    // Do a bounding box check so we're not doing unnecessary PointLocator lookups
    if (bbox.contains_point(p))
    {
      // First find the element the hit lands in
      const Elem * elem = getLocalElemContainingPoint(p, true);

      if (elem)
      {
        // We have to pass a vector of points into reinitElemPhys
        point_vec[0] = p;

        _subproblem.setCurrentSubdomainID(elem, 0);
        _subproblem.reinitElemPhys(elem, point_vec, 0); // Zero is for tid

        _values_positive_level_set_side[i] = (dynamic_cast<MooseVariable *>(_coupled_moose_vars[0]))
                                                 ->sln()[0]; // The zero is for the "qp"
        _grad_values_positive_level_set_side[i] =
            ((dynamic_cast<MooseVariable *>(_coupled_moose_vars[0]))->gradSln())[0];

        _found_points[i] = true;
      }

      const Elem * elem2 = getLocalElemContainingPoint(p, false);
      if (elem2)
      {
        // We have to pass a vector of points into reinitElemPhys
        point_vec[0] = p;

        _subproblem.setCurrentSubdomainID(elem2, 0);
        _subproblem.reinitElemPhys(elem2, point_vec, 0); // Zero is for tid

        _values_negative_level_set_side[i] = (dynamic_cast<MooseVariable *>(_coupled_moose_vars[0]))
                                                 ->sln()[0]; // The zero is for the "qp"
        _grad_values_negative_level_set_side[i] =
            ((dynamic_cast<MooseVariable *>(_coupled_moose_vars[0]))->gradSln())[0];

        _found_points[i] = true;
      }
    }
  }
}

void
PointValueAtXFEMInterface::finalize()
{
  // Save off for speed
  const auto pid = processor_id();

  /*
   * Figure out which processor is actually going "claim" each point.
   * If multiple processors found the point and computed values what happens is that
   * maxloc will give us the smallest PID in max_id
   */
  std::vector<unsigned int> max_id(_found_points.size());

  _communicator.maxloc(_found_points, max_id);

  for (auto i = beginIndex(max_id); i < max_id.size(); ++i)
  {
    // Only do this check on the proc zero because it's the same on every processor
    // _found_points should contain all 1's at this point (ie every point was found by a proc)
    if (pid == 0 && !_found_points[i])
      mooseError("In ", name(), ", sample point not found: ", _points[i]);

    // if (max_id[i] == pid)
    //   SamplerBase::addSample(_points[i], _ids[i], _point_values[i]);
  }

  // SamplerBase::finalize();
}

const Elem *
PointValueAtXFEMInterface::getLocalElemContainingPoint(const Point & p, bool positive_level_set)
{
  const Elem * elem1 = (*_pl)(p);
  // const Elem * elem = elem1;

  const Node * node = elem1->node_ptr(0);

  dof_id_type ls_dof_id = node->dof_number(_system.number(), _level_set_var_number, 0);
  Number ls_node_value = (*_solution)(ls_dof_id);

  bool positive = false;

  if (_xfem->isPointInsidePhysicalDomain(elem1, *node))
  {
    if (ls_node_value > 0.0)
      positive = true;
  }
  else
  {
    if (ls_node_value < 0.0)
      positive = true;
  }

  const Elem * elem2;

  for (auto & pair : *_elem_pairs)
  {
    if (pair.first == elem1)
      elem2 = pair.second;
    else if (pair.second == elem1)
      elem2 = pair.first;
  }

  if (((positive && positive_level_set) || (!positive && !positive_level_set)) &&
      (elem1 && elem1->processor_id() == processor_id()))
    return elem1;
  else if (((!positive && positive_level_set) || (positive && !positive_level_set)) &&
           (elem2 && elem2->processor_id() == processor_id()))
    return elem2;

  return nullptr;
}
