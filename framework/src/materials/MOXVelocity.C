/*************************************************/
/*           DO NOT MODIFY THIS HEADER           */
/*                                               */
/*                     BISON                     */
/*                                               */
/*    (c) 2015 Battelle Energy Alliance, LLC     */
/*            ALL RIGHTS RESERVED                */
/*                                               */
/*   Prepared by Battelle Energy Alliance, LLC   */
/*     Under Contract No. DE-AC07-05ID14517      */
/*     With the U. S. Department of Energy       */
/*                                               */
/*     See COPYRIGHT for full restrictions       */
/*                                               */
/*************************************************/

#include "MOXVelocity.h"
#include <iostream>
#include <string>


template<>
InputParameters validParams<MOXVelocity>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredCoupledVar("temp", "Coupled Temperature");
  params.addParam<Real>("limit", 1.0, "limit for pore velocity in m/s");
  params.addParam<Real>("pore_velocity", "calculated in materials property");

  return params;
}


MOXVelocity::MOXVelocity(const  InputParameters & parameters)
  :Material(parameters),
  _temp(coupledValue("temp")),
  _grad_temp(coupledGradient("temp")),
  _limit(getParam<Real>("limit")),
  _pore_velocity(declareProperty<Real>("pore_velocity"))
{}

MOXVelocity::~MOXVelocity()
{
}


void MOXVelocity::initQpStatefulProperties()
{
  _pore_velocity[_qp] = 0.0;
}


void MOXVelocity::computeQpProperties()
{
  _pore_velocity[_qp] = 1536.0*std::exp(-71763.0/_temp[_qp]) / (3.31e-24 * std::pow(_temp[_qp],1.5)) * 1.0;
  //    Moose::out << "pore speed = " << _pore_velocity[_qp];

  if (_pore_velocity[_qp] > _limit)
    _pore_velocity[_qp] = _limit;
}

