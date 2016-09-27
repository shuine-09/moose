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

#ifndef FISSIONRATEAUXLWR_H
#define FISSIONRATEAUXLWR_H

#include "AuxKernel.h"

class FissionRateAuxLWR : public AuxKernel
{
public:

  FissionRateAuxLWR(const InputParameters & parameters);

  virtual ~FissionRateAuxLWR() {}

protected:
  virtual Real computeValue();

  const Real _value;

  Function * const _function1;
  Function * const _function2;
  Function * const _function3;
  const Real _pellet_diameter;
  const Real _pellet_inner_diameter;
  const Real _energy_per_fission;
  const Real _fuel_volume_ratio;

  const Real _conversion_factor;
};

template<>
InputParameters validParams<FissionRateAuxLWR>();

#endif //FISSIONRATEAUXLWR_H
