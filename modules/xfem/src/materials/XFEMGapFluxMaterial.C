/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/
#include "XFEMGapFluxMaterial.h"

template<>
InputParameters validParams<XFEMGapFluxMaterial>()
{
  InputParameters params = validParams<Material>();
  params.addParam<Real>("heat_transfer_coef", 1.0, "heat transfer coefficient.");
  params.addParam<Real>("gap_tolerance", 1.0e-3, "gap tolerance.");
  params.addRequiredCoupledVar("temp", "Coupled Temperatrue");
  params.addRequiredCoupledVar("disp_x", "Displacement x");
  params.addRequiredCoupledVar("disp_y", "Displacement y");
  return params;
}

XFEMGapFluxMaterial::XFEMGapFluxMaterial(const InputParameters & parameters) :
    Material(parameters),
    _mat_prop(declareProperty<Real>("heatflux")),
    //_mat_prop_old(declarePropertyOld<Real>("heatflux")),
    _heat_transfer_coef(getParam<Real>("heat_transfer_coef")),
    _gap_tolerance(getParam<Real>("gap_tolerance")),
    _temp(coupledValue("temp")),
    _temp_neighbor(coupledNeighborValue("temp")),
    _disp_x(coupledValue("disp_x")),
    _disp_x_neighbor(coupledNeighborValue("disp_x")),
    _disp_y(coupledValue("disp_y")),
    _disp_y_neighbor(coupledNeighborValue("disp_y"))
{
}

void
XFEMGapFluxMaterial::computeQpProperties()
{
//  if (_temp_neighbor.size() == _temp.size())
  {
    Real gap = std::abs(_disp_y[_qp] - _disp_y_neighbor[_qp]);

    if (gap < _gap_tolerance)
      gap = _gap_tolerance;

    _mat_prop[_qp] = _heat_transfer_coef / gap * (_temp[_qp] - _temp_neighbor[_qp]);
    std::cout << "temp = " << _temp[_qp] << ", temp_neighbor = " << _temp_neighbor[_qp] << std::endl;  
    std::cout << "_mat_prop[_qp] = " << _mat_prop[_qp]  << std::endl;
  }
//  else
//    return;
  //  _mat_prop[_qp] = _heat_transfer_coef * _temp[_qp];
}
