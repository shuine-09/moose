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

#ifndef XFEMLEVELSET_H
#define XFEMLEVELSET_H

#include "AuxKernel.h"

class Function;

//Forward Declarations
class XFEMLevelSet;

template<>
InputParameters validParams<XFEMLevelSet>();

class XFEMLevelSet : public AuxKernel
{
public:

  XFEMLevelSet(const InputParameters & parameters);

  virtual ~XFEMLevelSet() {}

protected:
  virtual Real computeValue();

  Function & _func;
};

#endif // XFEMLEVELSET_H
