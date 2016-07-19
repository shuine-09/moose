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

#ifndef XFEMDIRICHLETBC_H
#define XFEMDIRICHLETBC_H

#include "NodalBC.h"

class XFEMDirichletBC;

template<>
InputParameters validParams<XFEMDirichletBC>();

/**
 * Boundary condition of a Dirichlet type
 *
 * Sets the value in the node
 */
class XFEMDirichletBC : public NodalBC
{
public:
  XFEMDirichletBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual bool shouldApply();
 
  /// The value for this BC
  const Real & _value;
};

#endif /* XFEMDIRICHLETBC_H */
