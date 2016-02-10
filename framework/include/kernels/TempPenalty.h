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

#ifndef TEMPPENALTY_H
#define TEMPPENALTY_H

#include "Kernel.h"

class TempPenalty;

template<>
InputParameters validParams<TempPenalty>();


class TempPenalty : public Kernel
{
public:
  TempPenalty(const InputParameters & parameters);
  virtual ~TempPenalty();

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

  Real _velocity;
  Real _diffusion;
  Real _cold_temp;
  Real _delta_temp;
  Real _cold_bdry;
  Real _hot_bdry;
  Real _alpha;

};


#endif /* TEMPPENALTY_H */
