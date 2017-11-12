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

#ifndef VESICLETOTALENERGYPOSTPROCESSOR_H
#define VESICLETOTALENERGYPOSTPROCESSOR_H

#include "ElementIntegralPostprocessor.h"
#include "MooseVariableInterface.h"

// Forward Declarations
class VesicleTotalEnergyPostprocessor;

template <>
InputParameters validParams<VesicleTotalEnergyPostprocessor>();

class VesicleTotalEnergyPostprocessor : public ElementIntegralPostprocessor,
                                        public MooseVariableInterface<Real>
{
public:
  VesicleTotalEnergyPostprocessor(const InputParameters & parameters);
  virtual void threadJoin(const UserObject & y);

protected:
  virtual Real computeQpIntegral();
  /// Holds the solution at current quadrature points
  const VariableValue & _u;
  /// Holds the solution gradient at the current quadrature points
  const VariableGradient & _grad_u;

  const VariableSecond & _second_u;

  Real _C;
  Real _epsilon;
  bool _rz;

  unsigned int _qp;
};

#endif
