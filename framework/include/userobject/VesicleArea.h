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

#ifndef VESICLEAREA_H
#define VESICLEAREA_H

#include "ElementUserObject.h"
#include "MooseVariableInterface.h"

// Forward Declarations
class VesicleArea;

template <>
InputParameters validParams<VesicleArea>();

class VesicleArea : public ElementUserObject, public MooseVariableInterface<Real>
{
public:
  VesicleArea(const InputParameters & parameters);

  virtual void initialize();
  virtual void execute();
  virtual void threadJoin(const UserObject & y);
  virtual void finalize();
  Real getValue() const;

protected:
  virtual Real computeQpIntegral();
  virtual Real computeIntegral();

  /// Holds the solution at current quadrature points
  const VariableValue & _u;
  /// Holds the solution gradient at the current quadrature points
  const VariableGradient & _grad_u;

  Real _epsilon;

  unsigned int _qp;

  Real _integral_value;
};

#endif
