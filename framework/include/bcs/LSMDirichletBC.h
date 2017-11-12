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

#ifndef LSMDIRICHLETBC_H
#define LSMDIRICHLETBC_H

#include "IntegratedBC.h"

//Forward Declarations
class LSMDirichletBC;

template<>
InputParameters validParams<LSMDirichletBC>();

class LSMDirichletBC : public IntegratedBC
{
public:

  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  LSMDirichletBC(const InputParameters & parameters);

  virtual ~LSMDirichletBC() {}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  const VariablePhiSecond & _second_phi;
  const VariableTestSecond & _second_test;
  const VariableSecond & _second_u;

private:
  Real _eps_m;
  Real _alpha;
};

#endif //LSMDIRICHLETBC_H
