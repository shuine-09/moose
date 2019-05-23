//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ConvectiveFluxFunction.h"
#include "libmesh/elem.h"
#include "Function.h"

registerMooseObject("HeatConductionApp", ConvectiveFluxFunction);

InputParameters
ConvectiveFluxFunction::validParams()
{
  InputParameters params = IntegratedBC::validParams();
  params.addRequiredParam<FunctionName>("T_infinity", "Function describing far-field temperature");
  params.addRequiredParam<FunctionName>("coefficient",
                                        "Function describing heat transfer coefficient");
  MooseEnum coef_func_type("TIME_AND_POSITION TEMPERATURE", "TIME_AND_POSITION");
  params.addParam<MooseEnum>(
      "coefficient_function_type",
      coef_func_type,
      "Type of function for heat transfer coefficient provided in 'coefficient' parameter");
  params.addClassDescription(
      "Determines boundary value by fluid heat transfer coefficient and far-field temperature");
  params.addParam<UserObjectName>("marker_uo", "Marker UserObject");

  return params;
}

ConvectiveFluxFunction::ConvectiveFluxFunction(const InputParameters & parameters)
  : IntegratedBC(parameters),
    _T_infinity(getFunction("T_infinity")),
    _coefficient(getFunction("coefficient")),
    _coef_func_type(getParam<MooseEnum>("coefficient_function_type").getEnum<CoefFuncType>()),
    _marker_uo(isParamValid("marker_uo") ? &getUserObjectByName<ActivatedElementsMarkerUO>(
                                               getParam<UserObjectName>("marker_uo"))
                                         : nullptr)
{
  if (_marker_uo)
    _marker_map = &(_marker_uo->getActivatedElementsMap());
  else
    _marker_map = nullptr;
}

Real
ConvectiveFluxFunction::computeQpResidual()
{
  bool apply = false;
  if (_marker_map)
  {
    dof_id_type elem_id = _current_elem->id();
    Real activate_elem = _marker_map->find(elem_id)->second;
    const Elem * neighbor_elem = _current_elem->neighbor_ptr(_current_side);
    if (neighbor_elem == nullptr && activate_elem > 0.5)
      apply = true;
    else if (neighbor_elem != nullptr && (activate_elem > 0.5) &&
             ((_marker_map->find(neighbor_elem->id())->second) < 0.5))
      apply = true;
  }

  Real coef;
  if (_coef_func_type == CoefFuncType::TIME_AND_POSITION)
    coef = _coefficient.value(_t, _q_point[_qp]);
  else
    coef = _coefficient.value(_u[_qp], Point());

  if (apply || _current_elem->subdomain_id() != 1)
    return _test[_i][_qp] * coef * (_u[_qp] - _T_infinity.value(_t, _q_point[_qp]));
  else
    return 0.0;
}

Real
ConvectiveFluxFunction::computeQpJacobian()
{
  bool apply = false;
  if (_marker_map)
  {
    dof_id_type elem_id = _current_elem->id();
    Real activate_elem = _marker_map->find(elem_id)->second;
    const Elem * neighbor_elem = _current_elem->neighbor_ptr(_current_side);
    if (neighbor_elem == nullptr && activate_elem > 0.5)
      apply = true;
    else if (neighbor_elem != nullptr && (activate_elem > 0.5) &&
             ((_marker_map->find(neighbor_elem->id())->second) < 0.5))
      apply = true;
  }

  Real coef;
  if (_coef_func_type == CoefFuncType::TIME_AND_POSITION)
    coef = _coefficient.value(_t, _q_point[_qp]);
  else
    coef = _coefficient.value(_u[_qp], Point());

  if (apply || _current_elem->subdomain_id() != 1)
    return _test[_i][_qp] * coef * _phi[_j][_qp];
  else
    return 0.0;
}
