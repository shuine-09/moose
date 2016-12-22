/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef MATERIALTENSORINTEGRALXFEM_H
#define MATERIALTENSORINTEGRALXFEM_H

#include "ElementIntegralPostprocessor.h"
#include "RankTwoTensor.h"

//Forward Declarations
class MaterialTensorIntegralXFEM;

class XFEM;

template<>
InputParameters validParams<MaterialTensorIntegralXFEM>();

/**
 * This postprocessor computes an element integral of a
 * component of a material tensor as specified by the user-supplied indices.
 */
class MaterialTensorIntegralXFEM: public ElementIntegralPostprocessor
{
public:
  MaterialTensorIntegralXFEM(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral();
  virtual Real computeIntegral();
  virtual Real getValue();

private:
  const MaterialProperty<RankTwoTensor> & _tensor;
  const unsigned int _i;
  const unsigned int _j;
  MooseSharedPointer<XFEM> _xfem;

};

#endif //MATERIALTENSORINTEGRALXFEM_H
