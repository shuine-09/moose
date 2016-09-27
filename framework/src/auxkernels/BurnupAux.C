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

#include "BurnupAux.h"

#include "Material.h"

template<>
InputParameters validParams<BurnupAux>()
{
  InputParameters params = validParams<AuxKernel>();

  params.addRequiredCoupledVar("fission_rate", "Coupled Fission Rate");
  params.addRequiredParam<Real>("density", "Initial fuel density");
  params.addParam<Real>("molecular_weight", 0.270, "The molecular weight");
  return params;
}

BurnupAux::BurnupAux(const InputParameters & parameters)
  :AuxKernel(parameters),
   _molecular_weight(getParam<Real>("molecular_weight")),
   _density(getParam<Real>("density")),
   _fission_rate(coupledValue("fission_rate"))
{}

Real
BurnupAux::computeValue()
{
  double avagadro_number = 6.023e23;

  double number_heavy_atoms = _density * avagadro_number/_molecular_weight;

  return _u_old[_qp] + (_fission_rate[_qp]*_dt)/number_heavy_atoms;
}
