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

#ifndef NEUTRONHEATSOURCE_H
#define NEUTRONHEATSOURCE_H

#include "FuelPinGeometry.h"
#include "Function.h"
#include "Kernel.h"

//Forward Declarations
class BurnupFunction;

class NeutronHeatSource : public Kernel
{
public:

  NeutronHeatSource(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

private:
  const Real _energy_per_fission;
  const std::string _units;
  const bool _has_fission_rate;
  const VariableValue & _fission_rate;

  const bool _has_max_fission_rate;
  const VariableValue & _max_fission_rate;
  const Real & _decay_heat_function;

  BurnupFunction * const _burnup_function;

  Function * const _rod_ave_lin_pow;
  Function * const _axial_profile;
  const FuelPinGeometry * const _pin_geometry;
  const bool _has_diameter;
  Real _outer_diameter;
  Real _inner_diameter;
  const bool _has_area;
  Real _area;
  const Real _fraction;

};

template<>
InputParameters validParams<NeutronHeatSource>();

#endif
