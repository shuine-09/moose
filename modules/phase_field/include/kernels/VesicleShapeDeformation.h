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

#ifndef VESICLESHAPEDEFORMATION_H
#define VESICLESHAPEDEFORMATION_H

#include "Kernel.h"

class VesicleShapeDeformation;

template <>
InputParameters validParams<VesicleShapeDeformation>();

class VesicleShapeDeformation : public Kernel
{
public:
  VesicleShapeDeformation(const InputParameters & parameters);

  virtual ~VesicleShapeDeformation();

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  Real _epsilon;
  Real _C;
  bool _rz;

  const VariablePhiSecond & _second_phi;
  const VariableTestSecond & _second_test;
  const VariableSecond & _second_u;
};

#endif /* VESICLESHAPEDEFORMATION_H */
