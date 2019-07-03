//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MaterialRealDenseVectorAux.h"

registerMooseObject("MooseApp", MaterialRealDenseVectorAux);

template <>
InputParameters
validParams<MaterialRealDenseVectorAux>()
{
  InputParameters params = validParams<MaterialAuxBase<>>();
  params.addParam<unsigned int>("index", 0, "The index to consider for this kernel");
  return params;
}

MaterialRealDenseVectorAux::MaterialRealDenseVectorAux(const InputParameters & parameters)
  : MaterialAuxBase<DenseVector<Real>>(parameters), _index(getParam<unsigned int>("index"))
{
}

Real
MaterialRealDenseVectorAux::getRealValue()
{
  return _prop[_qp](_index);
}
