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

#ifndef HEATTESTAUX_H
#define HEATTESTAUX_H

#include "AuxKernel.h"

//Forward Declarations
class HeatTestAux;

template<>
InputParameters validParams<HeatTestAux>();

/**
 * Function auxiliary value
 */
class HeatTestAux : public AuxKernel
{
public:
  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  HeatTestAux(const InputParameters & parameters);

  virtual ~HeatTestAux() {}

protected:
  virtual Real computeValue();

  Real _velocity;
  Real _diffusion;
  Real _cold_temp;
  Real _delta_temp;
  Real _cold_bdry;
  Real _hot_bdry;

  Real step_func(Real x);
};

#endif // HEATTESTAUX_H
