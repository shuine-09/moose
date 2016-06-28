/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "RankTwoStdVectorAux.h"
#include "RankTwoScalarTools.h"

template<>
InputParameters validParams<RankTwoStdVectorAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("Access a component of a RankTwoTensor stored in a std vector");
  params.addRequiredParam<MaterialPropertyName>("rank_two_tensor_vector", "The rank two material tensor name");
  params.addParam<unsigned int>("index", 0, "The index to consider for this kernel");
  params.addRequiredRangeCheckedParam<unsigned int>("index_i", "index_i >= 0 & index_i <= 2", "The index i of ij for the tensor to output (0, 1, 2)");
  params.addRequiredRangeCheckedParam<unsigned int>("index_j", "index_j >= 0 & index_j <= 2", "The index j of ij for the tensor to output (0, 1, 2)");
  return params;
}

RankTwoStdVectorAux::RankTwoStdVectorAux(const InputParameters & parameters) :
    AuxKernel(parameters),
    _tensors(getMaterialProperty<std::vector<RankTwoTensor> >("rank_two_tensor_vector")),
    _index(getParam<unsigned int>("index")),
    _i(getParam<unsigned int>("index_i")),
    _j(getParam<unsigned int>("index_j"))
{
}

Real
RankTwoStdVectorAux::computeValue()
{
  mooseAssert(_tensors[_qp].size() > _index, "RankTwoStdVectorAux: You chose to extract component " << _index << " but your RankTwoTensor vector only has size " << _tensors[_qp].size());
  return RankTwoScalarTools::component(_tensors[_qp][_index], _i, _j);
}
