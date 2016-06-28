/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef RANKTWOSTDVECTORAUX_H
#define RANKTWOSTDVECTORAUX_H

#include "AuxKernel.h"
#include "RankTwoTensor.h"

class RankTwoStdVectorAux;

/**
 * RankTwoStdVectorAux is designed to take the data in the RankTwoTensor material
 * property vector, for example stress or strain, and output the value for the
 * supplied indices.
 */

template<>
InputParameters validParams<RankTwoStdVectorAux>();

class RankTwoStdVectorAux : public AuxKernel
{
public:
  RankTwoStdVectorAux(const InputParameters & parameters);

protected:
  virtual Real computeValue();

private:
  const MaterialProperty<std::vector<RankTwoTensor> > & _tensors;
  
  /// index of the vecor element
  unsigned int _index;
  
  const unsigned int _i;
  const unsigned int _j;
};

#endif //RANKTWOSTDVECTORAUX_H
