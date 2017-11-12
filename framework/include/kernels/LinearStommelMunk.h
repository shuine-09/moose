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

#ifndef LINEARSTOMMELMUNK_H
#define LINEARSTOMMELMUNK_H

#include "Kernel.h"

class LinearStommelMunk;

template<>
InputParameters validParams<LinearStommelMunk>();


class LinearStommelMunk : public Kernel
{
public:
  LinearStommelMunk(const InputParameters & parameters);
  virtual ~LinearStommelMunk();

protected:
  const VariablePhiSecond & _second_phi;
  const VariableTestSecond & _second_test;
  const VariableSecond & _second_u;

  Real _eps_s;
  Real _eps_m;

  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
};


#endif /* LINEARSTOMMELMUNK_H */
