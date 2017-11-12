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

#ifndef VESICLEAREAPOSTPROCESSOR_H
#define VESICLEAREAPOSTPROCESSOR_H

#include "ElementIntegralPostprocessor.h"
#include "MooseVariableInterface.h"

// Forward Declarations
class VesicleAreaPostprocessor;

template <>
InputParameters validParams<VesicleAreaPostprocessor>();

class VesicleAreaPostprocessor : public ElementIntegralPostprocessor,
                                 public MooseVariableInterface<Real>
{
public:
  VesicleAreaPostprocessor(const InputParameters & parameters);
  virtual void threadJoin(const UserObject & y);

protected:
  virtual Real computeQpIntegral();
  /// Holds the solution at current quadrature points
  const VariableValue & _u;
  /// Holds the solution gradient at the current quadrature points
  const VariableGradient & _grad_u;

  Real _epsilon;

  unsigned int _qp;
};

#endif
