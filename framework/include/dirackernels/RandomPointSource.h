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

#ifndef RANDOMPOINTSOURCE_H
#define RANDOMPOINTSOURCE_H

// Moose Includes
#include "DiracKernel.h"

//Forward Declarations
class RandomPointSource;

template<>
InputParameters validParams<RandomPointSource>();

/**
 * Insert randomly placed point impulses into the simulatiom domain.
 * \todo Should probably use the new RandomInterface
 */
class RandomPointSource : public DiracKernel
{
public:
  RandomPointSource(const InputParameters & parameters);

  virtual void addPoints();
  virtual Real computeQpResidual();

protected:
  const unsigned int _mesh_dimension;

  Real _value;
  unsigned int _num_srcs;
  Real _dt_one;
  Point _p;

  Real _x1, _y1, _z1;
  Real _x2, _y2, _z2;

  Point _bottom_left, _top_right, _range;
};

#endif //RANDOMPOINTSOURCE_H
