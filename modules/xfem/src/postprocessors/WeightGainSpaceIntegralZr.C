//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "WeightGainSpaceIntegralZr.h"

registerMooseObject("MooseApp", WeightGainSpaceIntegralZr);

defineLegacyParams(WeightGainSpaceIntegralZr);

InputParameters
WeightGainSpaceIntegralZr::validParams()
{
  InputParameters params = GeneralPostprocessor::validParams();
  params.addParam<Real>("temperature",1473.15,"Temperature of the cladding (K)");
  params.addRequiredParam<PostprocessorName>("concentration_integral", "The name of the oxygen concentration integral postprocessor");
  params.addParam<Real>("ymax",4,"The mesh y dimension [um]");
  params.addRequiredParam<PostprocessorName>("oxide_thickness", "The name of the oxide thickness postprocessor");
  params.addRequiredParam<PostprocessorName>("alpha_thickness", "The name of the alpha layer thickness postprocessor");
  return params;
}

WeightGainSpaceIntegralZr::WeightGainSpaceIntegralZr(const InputParameters & parameters)
  : GeneralPostprocessor(parameters),
    _temperature(getParam<Real>("temperature")),
    _wg(0),
    _C_integral(getPostprocessorValue("concentration_integral")),
    _ymax(getParam<Real>("ymax")),
    _delta(getPostprocessorValue("oxide_thickness")),
    _d_alpha(getPostprocessorValue("alpha_thickness"))
{
}

void
WeightGainSpaceIntegralZr::initialize()
{
}

void
WeightGainSpaceIntegralZr::execute()
{
  const Real Mo = 15.99;    // Molar mass of oxygen
  const Real Na(6.022e23);
  const Real Mzr(91.22);
  const Real rho_Zr(6.56);
  const Real Czr = 1e6 * rho_Zr * Na / Mzr ; // Concentration of Zr atoms in the metal [/m³]
  const Real rho_ZrO2(5.68);
  const Real Czro2 = 1e6 * rho_ZrO2 * Na / (Mzr+2*Mo) ; // Concentration of Zr atoms in the oxide[/m³]
  const Real Zr_PBR(1.55);

  // First, convert the integrated concentration in /m² (divide by ymax to cancel integration over y)
  Real full_weak_integral = _C_integral / _ymax * Czr * 1e-6; ///!!!units!!!

  // Then, remove oxide part of the integral.
  // for that, need various interfaces concentrations
  const Real x_ox = 0.66666666667;
  const Real C_ox = 3 * Czro2 * x_ox;

  const Real x_ox_a = 0.667118123 - 1.10606e-5 * _temperature;
  const Real C_ox_a = 3 * Czro2 * x_ox_a;

  Real x_a_ox;
  if (_temperature > 473.15 && _temperature < 1478.15)
  {
    x_a_ox = (28.6 + exp(-6748/_temperature + 4.748)) * 1e-2;
  }
  else if (_temperature > 1478.15 && _temperature < 1798.15)
  {
    x_a_ox = (28.6 + exp(-6301/_temperature + 4.460)) * 1e-2;
  }
  else if (_temperature > 1798.15 && _temperature < 2338.15)
  {
    x_a_ox = (28.6 + exp(-7012/_temperature + 8.434 - 3.521e-3 * _temperature )) * 1e-2;
  }
  else
  {
    x_a_ox = 28.6 *1e-2;
  }
  const Real C_a_ox = Czr * x_a_ox/(1-x_a_ox);

  Real x_b_a = (9.59e-3 * (_temperature - 1136) + 4.72e-6 * pow(_temperature - 1136,2) - 4.35e-9 * pow(_temperature - 1136,3)) * 1e-2;
  Real x_a_b = (45.86e-3 * (_temperature - 1136) - 44.77e-6 * pow(_temperature - 1136,2) + 17.40e-9 * pow(_temperature - 1136,3)) * 1e-2; //the original one, not the weak equivalent
  const Real C_a_b = Czr * x_a_b / (1 - x_a_b);
  const Real C_b_a = Czr * x_b_a / (1 - x_b_a);

  // can now get the oxide part of the weak integral [/m²]:
  Real oxide_weak_integral = 1e-6 * _delta/Zr_PBR * (C_a_ox - (C_a_b-C_b_a) + 0.5*(C_ox-C_ox_a));

  // And from that the metal part only of the weak integral [/m²]:
  Real metal_weak_integral = full_weak_integral - oxide_weak_integral;

  // The strong part of the integral in the metal, that was not accounted for during integration [/m²]:
  Real metal_strong_missing = 1e-6 * _d_alpha * (C_a_b-C_b_a);

  // The integral in the oxide (linear concentration) [/m²]:
  Real oxide_strong_integral = 1e-6 * _delta * 0.5*(C_ox+C_ox_a);

  // we also have to remove the initial oxygen natively in the metal [/m²]:
  const Real x_nat = 0.0074;
  const Real C_nat = Czr * x_nat/(1-x_nat);
  Real metal_native_integral = 1e-6 * 600 * C_nat;

  // Thus, the integral of the oxygen concentration that has came into the sample during corrosion is [/m²]:
  Real total_integral = metal_weak_integral + metal_strong_missing + oxide_strong_integral - metal_native_integral;

//  std::cout << "full_weak_integral : " << full_weak_integral << std::endl;
//  std::cout << "oxide_weak_integral : " << oxide_weak_integral << std::endl;
//  std::cout << "metal_weak_integral : " << metal_weak_integral << std::endl;
//  std::cout << "metal_strong_missing : " << metal_strong_missing << std::endl;
//  std::cout << "oxide_strong_integral : " << oxide_strong_integral << std::endl;
//  std::cout << "metal_native_integral : " << metal_native_integral << std::endl;
//  std::cout << "total_integral : " << total_integral << std::endl;

  // Finally, we can get the weight gain resulting from space integration of the oxygen concentration [mg/cm²]:
 _wg = 0.1 * Mo / Na * total_integral;
}

Real
WeightGainSpaceIntegralZr::getValue()
{
    return _wg;
}
