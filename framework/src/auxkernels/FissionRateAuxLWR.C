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

/*FissionRateAuxLWR is the auxkernel used to calculate the fission rate.  It is specific to light water reactors.
 * PiecewiseBilinear is used to calculate the axial power profile, which is a funciton of axial position (typically the y-coordinate of the node)
 * and time, function2 here.  The fission rate is calculated here by taking the rod average linear power (function1 here)
 * multiplying it by the axial power profile, dividing by the cross sectional area of the pellet, dividing by the energy per fission,
 * and multipling buy a volume reduction term,
 * f_volume_recduction, which accounts for the dishing of the pellet (i.e.deviation from right circular cylinder geometry.  The fission rate
 * is used to calculate the source term in NeutronHeatSource and the burnup in BurnupAux.
 *
 * The use of the variable called value is continued from AuxKernel.  It's use is continued for debugging purposes and to simplify the development of classes
 * that inherit from Auxkernel.  The variable value can be specified in the input file, for which case it would act as a scale factor.
 *
 * */

#include "FissionRateAuxLWR.h"
#include "Function.h"

template<>
InputParameters validParams<FissionRateAuxLWR>()
{
  InputParameters params = validParams<AuxKernel>();
  params.set<Real>("value")=1.0;
  params.addParam<FunctionName>("rod_ave_lin_pow", "", "The rod power history function.");
  params.addParam<FunctionName>("axial_power_profile", "", "The axial power profile function.");
  params.addParam<FunctionName>("radial_power_profile", "", "The radial power profile function.");
  params.addRequiredParam<Real>("pellet_diameter", "Fuel pellet diameter in m");
  params.addParam<Real>("pellet_inner_diameter", 0, "Pellet inner diameter in m for an annular pellet");
  params.addParam<Real>("energy_per_fission", 3.28451e-11, "energy per fission in J/fission");
  params.addParam<Real>("fuel_volume_ratio", 1, "Reduction factor for deviation from right circular cylinder fuel.  The ratio of actual volume to right circular cylinder volume.");
//  params.set<bool>("use_displaced_mesh") = true;
  return params;
}

FissionRateAuxLWR::FissionRateAuxLWR(const InputParameters & parameters)
  :AuxKernel(parameters),
   _value(getParam<Real>("value")),
   _function1( getParam<FunctionName>("rod_ave_lin_pow") != "" ? &getFunction("rod_ave_lin_pow") : NULL ),
   _function2( getParam<FunctionName>("axial_power_profile") != "" ? &getFunction("axial_power_profile") : NULL ),
   _function3( getParam<FunctionName>("radial_power_profile") != "" ? &getFunction("radial_power_profile") : NULL ),
   _pellet_diameter(parameters.get<Real>("pellet_diameter")),
   _pellet_inner_diameter(parameters.get<Real>("pellet_inner_diameter")),
   _energy_per_fission(parameters.get<Real>("energy_per_fission")),
   _fuel_volume_ratio(parameters.get<Real>("fuel_volume_ratio")),
   _conversion_factor( 4.0 / M_PI / (_pellet_diameter*_pellet_diameter - _pellet_inner_diameter*_pellet_inner_diameter) / _energy_per_fission / _fuel_volume_ratio )
{
  if (_fuel_volume_ratio > 1)
  {
    mooseError("FissionRateAuxLWR: fuel_volume_ratio must be <= 1");
  }
}

Real
FissionRateAuxLWR::computeValue()
{
  Real value( _value );
  if ( _function1 )
  {
    if ( isNodal() )
    {
      mooseAssert( _current_node, "ERROR:  FissionRateAuxLWR reports nodal but no node defined." );
      const Node & node = *_current_node;
      value *= _function1->value(_t, node);
      if (_function2)
      {
        value *= _function2->value(_t, node);
      }
      if (_function3)
     {
        value *= _function3->value(_t, node);
      }
    }
   else
    {
      value *= _function1->value(_t, _q_point[_qp]);
      if (_function2)
      {
        value *= _function2->value(_t, _q_point[_qp]);
      }
      if (_function3)
      {
        value *= _function3->value(_t, _q_point[_qp]);
      }
    }
    value *= _conversion_factor;
  }
  if (value < 0)
  {
    mooseError("Negative fission rate in FissionRateAuxLWR");
  }
  return value;
}
