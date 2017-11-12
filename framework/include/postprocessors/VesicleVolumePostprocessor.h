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

#ifndef VESICLEVOLUMEPOSTPROCESSOR_H
#define VESICLEVOLUMEPOSTPROCESSOR_H

#include "ElementIntegralPostprocessor.h"
#include "MooseVariableInterface.h"

// Forward Declarations
class VesicleVolumePostprocessor;

template <>
InputParameters validParams<VesicleVolumePostprocessor>();

class VesicleVolumePostprocessor : public ElementIntegralPostprocessor,
                                   public MooseVariableInterface<Real>
{
public:
  VesicleVolumePostprocessor(const InputParameters & parameters);
  virtual void threadJoin(const UserObject & y);

protected:
  virtual Real computeQpIntegral();
  /// Holds the solution at current quadrature points
  const VariableValue & _u;
  /// Holds the solution gradient at the current quadrature points
  const VariableGradient & _grad_u;

  unsigned int _qp;
};

#endif
