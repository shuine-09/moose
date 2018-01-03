/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "MaximumNormalSeparation.h"
#include "Assembly.h"

template <>
InputParameters
validParams<MaximumNormalSeparation>()
{
  InputParameters params = validParams<Material>();
  params.addClassDescription("");
  params.addRequiredCoupledVar("disp_x", "Name of the variable to couple");
  params.addRequiredCoupledVar("disp_y", "Name of the variable to couple");
  params.addParam<std::string>("base_name",
                               "Optional parameter that allows the user to define "
                               "multiple mechanics material systems on the same block");
  return params;
}

MaximumNormalSeparation::MaximumNormalSeparation(const InputParameters & parameters)
  : Material(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _max_normal_separation(declareProperty<Real>(_base_name + "max_normal_separation")),
    _max_normal_separation_old(getMaterialPropertyOld<Real>(_base_name + "max_normal_separation")),
    _disp_x(coupledValue("disp_x")),
    _disp_x_neighbor(coupledNeighborValue("disp_x")),
    _disp_y(coupledValue("disp_y")),
    _disp_y_neighbor(coupledNeighborValue("disp_y")),
    _normals(_assembly.normals())
{
}

void
MaximumNormalSeparation::computeQpProperties()
{
  _max_normal_separation[_qp] = _disp_y[_qp] - _disp_y_neighbor[_qp];
}

void
MaximumNormalSeparation::initQpStatefulProperties()
{
  _max_normal_separation[_qp] = 0.0;
}
