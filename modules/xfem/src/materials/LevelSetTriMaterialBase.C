//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LevelSetTriMaterialBase.h"
#include "AuxiliarySystem.h"
#include "MooseVariable.h"
#include "XFEM.h"

InputParameters
LevelSetTriMaterialBase::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Compute a material property for tri-materials (consisting of three "
                             "different materials) defined by a level set function.");
  params.addRequiredParam<VariableName>(
      "ls_var_1", "The name of level set variable used to represent the left interface");
  params.addRequiredParam<VariableName>(
      "ls_var_2", "The name of level set variable used to represent the right interface");
  params.addRequiredParam<std::string>("levelset_pos_pos_base",
                                       "Base name for the material in level set positive-positive region.");
  params.addRequiredParam<std::string>("levelset_pos_neg_base",
                                       "Base name for the material in level set positive-negative region.");
  params.addRequiredParam<std::string>("levelset_neg_neg_base",
                                       "Base name for the material in level set negative-negative region.");
  params.addParam<std::string>("base_name",
                               "Base name for the computed material property (optional)");
  params.addRequiredParam<std::string>("prop_name", "Name for the computed material property.");
  return params;
}

LevelSetTriMaterialBase::LevelSetTriMaterialBase(const InputParameters & parameters)
  : Material(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _prop_name(getParam<std::string>("prop_name")),
    _ls_var_1_number(_subproblem
                              .getVariable(_tid,
                                           parameters.get<VariableName>("ls_var_1"),
                                           Moose::VarKindType::VAR_ANY,
                                           Moose::VarFieldType::VAR_FIELD_STANDARD)
                              .number()),
    _ls_var_2_number(_subproblem
                              .getVariable(_tid,
                                           parameters.get<VariableName>("ls_var_2"),
                                           Moose::VarKindType::VAR_ANY,
                                           Moose::VarFieldType::VAR_FIELD_STANDARD)
                              .number()),
    _system1(_subproblem.getSystem(getParam<VariableName>("ls_var_1"))),
    _solution1(_system1.current_local_solution.get()),
    _system2(_subproblem.getSystem(getParam<VariableName>("ls_var_2"))),
    _solution2(_system2.current_local_solution.get()),
    _use_neg_neg_property(false),
    _use_pos_neg_property(false)
{
  FEProblemBase * fe_problem = dynamic_cast<FEProblemBase *>(&_subproblem);

  if (fe_problem == NULL)
    mooseError("Problem casting _subproblem to FEProblemBase in XFEMMaterialStateMarkerBase");

  _xfem = MooseSharedNamespace::dynamic_pointer_cast<XFEM>(fe_problem->getXFEM());
}

void
LevelSetTriMaterialBase::computeProperties()
{
  const Node * node = _current_elem->node_ptr(0);

  dof_id_type ls1_dof_id = node->dof_number(_system1.number(), _ls_var_1_number, 0);
  Number ls1_node_value = (*_solution1)(ls1_dof_id);

  dof_id_type ls2_dof_id = node->dof_number(_system2.number(), _ls_var_2_number, 0);
  Number ls2_node_value = (*_solution2)(ls2_dof_id);

  _use_neg_neg_property = false;
  _use_pos_neg_property = false;

  if (_xfem->isPointInsidePhysicalDomain(_current_elem, *node))
  {
    if (ls1_node_value < 0.0)
      _use_neg_neg_property = true;
    else
    {
      if (ls2_node_value < 0.0)
      _use_pos_neg_property = true;
    }
  }
  else //nonphysical nodes
  {
    if ((ls1_node_value < 0.0) || (ls2_node_value > 0.0))
      _use_pos_neg_property = true;
    else
    {
      if (abs(ls1_node_value) < abs(ls2_node_value)) //right side of the left interface
        _use_neg_neg_property = true;
    }
  }

  Material::computeProperties();
}

void
LevelSetTriMaterialBase::computeQpProperties()
{
  if (_use_neg_neg_property)
    assignQpPropertiesForLevelSetNegNeg();
  else if (_use_pos_neg_property)
    assignQpPropertiesForLevelSetPosNeg();
  else
    assignQpPropertiesForLevelSetPosPos();
}
