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
#include "StatefulHeatFluxMaterial.h"

template<>
InputParameters validParams<StatefulHeatFluxMaterial>()
{
  InputParameters params = validParams<Material>();
  params.addParam<Real>("heat_transfer_coef", 1.0, "heat transfer coefficient.");
  params.addRequiredCoupledVar("temp", "Coupled Temperatrue");
  return params;
}

StatefulHeatFluxMaterial::StatefulHeatFluxMaterial(const InputParameters & parameters) :
    Material(parameters),
    _mat_prop(declareProperty<Real>("heatflux")),
    _mat_prop_old(declarePropertyOld<Real>("heatflux")),
    _heat_transfer_coef(getParam<Real>("heat_transfer_coef")),
    _temp(coupledValue("temp")),
    _temp_neighbor(coupledNeighborValue("temp"))
{
}

void
StatefulHeatFluxMaterial::initQpStatefulProperties()
{
  _mat_prop[_qp] = 1.0;
}


void
StatefulHeatFluxMaterial::computeQpProperties()
{
  if (_temp_neighbor.size() == _temp.size())
  {
    //_mat_prop[_qp] = _heat_transfer_coef * (_temp[_qp] - _temp_neighbor[_qp]);
    _mat_prop[_qp] = _mat_prop_old[_qp]  + 1.0;
    std::cout << "temp = " << _temp[_qp] << ", temp_neighbor = " << _temp_neighbor[_qp] << std::endl;  
  }
  else 
    _mat_prop[_qp] = -1.0;
}
