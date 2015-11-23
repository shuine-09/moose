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

#ifndef XFEMWEIBULLAUX_H
#define XFEMWEIBULLAUX_H

#include "AuxKernel.h"

class XFEM;

/**
 * Coupled auxiliary value
 */
class XFEMWeibullAux : public AuxKernel
{
public:

  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  XFEMWeibullAux(const InputParameters & parameters);

  virtual ~XFEMWeibullAux() {}

protected:
  virtual Real computeValue();

private:
  XFEM *_xfem;
  int _weibull_modulus;
  Real _specimen_volume;
  Real _specimen_material_property;
};

template<>
InputParameters validParams<XFEMWeibullAux>();

#endif //XFEMWEIBULLAUX_H
