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

#ifndef TWOSIDEDDIFFUSION_H
#define TWOSIDEDDIFFUSION_H

#include "Kernel.h"

class XFEM;

class TwoSidedDiffusion;

template<>
InputParameters validParams<TwoSidedDiffusion>();

class TwoSidedDiffusion : public Kernel
{
public:
  TwoSidedDiffusion(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;

  bool isPlusDomain();

  MooseSharedPointer<XFEM> _xfem;

  Real _diffusivity_ls_plus;
  Real _diffusivity_ls_minus;

  const VariableValue & _ls;
};


#endif /* TWOSIDEDDIFFUSION_H */
