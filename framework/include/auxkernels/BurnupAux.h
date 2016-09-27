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

#ifndef BURNUPAUX_H
#define BURNUPAUX_H

#include "AuxKernel.h"

/**
 * Coupled auxiliary value
 */
class BurnupAux : public AuxKernel
{
public:

  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  BurnupAux(const InputParameters & parameters);

  virtual ~BurnupAux() {}

protected:
  virtual Real computeValue();

private:
  const Real _molecular_weight;
  const Real _density;
  const VariableValue & _fission_rate;
};

template<>
InputParameters validParams<BurnupAux>();

#endif //BURNUPAUX_H
