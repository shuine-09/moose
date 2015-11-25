/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
//  Crystal plasticity userobject system  base class.
//
#include "CrystalPlasticityUOBase.h"

template<>
InputParameters validParams<CrystalPlasticityUOBase>()
{
  InputParameters params = validParams<ElementUserObject>();
  params.addRequiredParam<std::string>("uo_name", "The userobject name.");
  params.addRequiredParam<unsigned int>("uo_data_size","The userobject's data size, i.e. the number of slip of systems.");
  params.addClassDescription("Crystal plasticity userobject system base class.  Override the virtual functions in your class");
  return params;
}

CrystalPlasticityUOBase::CrystalPlasticityUOBase(const InputParameters & parameters) :
  ElementUserObject(parameters),
  _uo_name(getParam<std::string>("uo_name")),
  _uo_data_size(getParam<unsigned int>("uo_data_size"))
{
}

std::string
CrystalPlasticityUOBase::crystalPlasticityUOName() const
{
  return _uo_name;
}

unsigned int
CrystalPlasticityUOBase::dataSize() const
{
  return _uo_data_size;
}
