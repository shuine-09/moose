/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef RANKTWOAUXXFEM_H
#define RANKTWOAUXXFEM_H

#include "AuxKernel.h"
#include "RankTwoTensor.h"

class XFEM;
class RankTwoAuxXFEM;

/**
 * RankTwoAuxXFEM is designed to take the data in the RankTwoTensor material
 * property, for example stress or strain, and output the value for the
 * supplied indices.
 */

template<>
InputParameters validParams<RankTwoAuxXFEM>();

class RankTwoAuxXFEM : public AuxKernel
{
public:
  RankTwoAuxXFEM(const InputParameters & parameters);

protected:
  virtual Real computeValue();
  virtual void compute(); 

private:
  const MaterialProperty<RankTwoTensor> & _tensor;
  const unsigned int _i;
  const unsigned int _j;
  MooseSharedPointer<XFEM> _xfem;

};

#endif //RANKTWOAUXXFEM_H
