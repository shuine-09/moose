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

#ifndef INTERFACEGLUEDCONTACT_H
#define INTERFACEGLUEDCONTACT_H

#include "InterfaceKernel.h"

// Forward Declarations
class InterfaceGluedContact;

template <>
InputParameters validParams<InterfaceGluedContact>();

/**
 * DG kernel for interfacing diffusion between two variables on adjacent blocks
 */
class InterfaceGluedContact : public InterfaceKernel
{
public:
  InterfaceGluedContact(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual(Moose::DGResidualType type);
  virtual Real computeQpJacobian(Moose::DGJacobianType type);

  const Real _alpha;
  const unsigned int _component;
  const bool _tangential;
};

#endif
