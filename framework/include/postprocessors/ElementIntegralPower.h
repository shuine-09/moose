/*************************************************/
/*           DO NOT MODIFY THIS HEADER           */
/*                                               */
/*                     BISON                     */
/*                                               */
/*    (c) 2015 Battelle Energy Alliance, LLC     */
/*            ALL RIGHTS RESERVED                */
/*                                               */
/*   Prepared by Battelle Energy Alliance, LLC   */
/*     Under Contract No. DE-AC07-05ID14517      */
/*     With the U. S. Department of Energy       */
/*                                               */
/*     See COPYRIGHT for full restrictions       */
/*************************************************/

#ifndef ELEMENTINTEGRALPOWER_H
#define ELEMENTINTEGRALPOWER_H

#include "ElementIntegralVariablePostprocessor.h"

//Forward Declarations
class BurnupFunction;

/**
 * This postprocessor computes an integral of the volumetric fission rate times the energy per fission
 * giving the power
 *
 */
class ElementIntegralPower: public ElementIntegralVariablePostprocessor
{
public:
  ElementIntegralPower(const InputParameters & parameters);
  virtual ~ElementIntegralPower() {}

protected:
  virtual Real computeQpIntegral();
  const Real _energy_per_fission;
  const bool _has_fission_rate;
  const VariableValue & _fission_rate;
  BurnupFunction * const _burnup_function;
};

template<>
InputParameters validParams<ElementIntegralPower>();

#endif
