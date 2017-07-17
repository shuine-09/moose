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

#ifndef DIFFUSIONSTRESS_H
#define DIFFUSIONSTRESS_H

#include "Kernel.h"

class DiffusionStress;

template <>
InputParameters validParams<DiffusionStress>();

/**
 * This kernel implements the Laplacian operator:
 * $\nabla u \cdot \nabla \phi_i$
 */
class DiffusionStress : public Kernel
{
public:
  DiffusionStress(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;

  const VariableValue & _pressure;
  const VariableGradient & _grad_pressure;
  Real _coef;
};

#endif /* DIFFUSIONSTRESS_H */
