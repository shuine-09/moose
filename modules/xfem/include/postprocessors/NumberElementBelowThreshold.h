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

#ifndef NUMBERELEMENTBELOWTHRESHOLD_H
#define NUMBERELEMENTBELOWTHRESHOLD_H

#include "ElementIntegralPostprocessor.h"

// Forward Declarations
class NumberElementBelowThreshold;

template <>
InputParameters validParams<NumberElementBelowThreshold>();

class NumberElementBelowThreshold : public ElementIntegralPostprocessor
{
public:
  NumberElementBelowThreshold(const InputParameters & parameters);

  virtual Real computeQpIntegral() override;
  virtual Real computeIntegral() override;

protected:
  const MaterialProperty<Real> & _weibull;
};

#endif // NUMBERELEMENTBELOWTHRESHOLD_H
