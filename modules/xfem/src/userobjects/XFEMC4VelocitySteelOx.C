//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "XFEMC4VelocitySteelOx.h"

registerMooseObject("XFEMApp", XFEMC4VelocitySteelOx);

template <>
InputParameters
validParams<XFEMC4VelocitySteelOx>()
{
  InputParameters params = validParams<XFEMMovingInterfaceVelocityBase>();
  params.addClassDescription(
      "Calculate the oxide/gas interface velocity for the 1 oxide layer (spinel) C4 model for steel corrosion.");
  return params;
}

XFEMC4VelocitySteelOx::XFEMC4VelocitySteelOx(const InputParameters & parameters)
  : XFEMMovingInterfaceVelocityBase(parameters)
{
}

Real
XFEMC4VelocitySteelOx::computeMovingInterfaceVelocity(unsigned int point_id) const
{
  //RealVectorValue grad_positive = _value_at_interface_uo->getGradientAtPositiveLevelSet()[point_id];
  //RealVectorValue grad_negative = _value_at_interface_uo->getGradientAtNegativeLevelSet()[point_id];

  Real xt = (_value_at_interface_uo->getPointCurrentLocation(point_id))(0);

  //Oxide thickness [nm]
  Real delta = std::abs(xt - 5000);
  std::cout << "delta : " << delta << std::endl;


  // Current implementation only supports the case that the interface is moving horizontally

  // Fundamental physics constants :
  const Real Na(6.022140857e23);
  const Real Kb(8.6173303e-5);

/*  // 21-2N steel composition (atomic fractions, [at%])

  const Real x_steel_Fe(0.6548);
  const Real x_steel_Cr(0.2100);
  const Real x_steel_Mn(0.0834);
  const Real x_steel_Ni(0.0196);
  const Real x_steel_C(0.0247);
  const Real x_steel_Mo(0.0028);
  const Real x_steel_Si(0.0048);
*/
  // Material properties

  // Atomic molar masses [g/mol]
  const Real M_Fe = 55.845;
  const Real M_Cr = 54.938044;
  const Real M_Mn = 51.9961;
  const Real M_Ni = 58.6934;
  const Real M_C = 12.0107;
  const Real M_Mo = 95.95;
  const Real M_Si = 28.0855;
//  const Real M_steel = x_steel_Fe * M_Fe + x_steel_Cr * M_Cr + x_steel_Mn * M_Mn + x_steel_Ni * M_Ni + x_steel_C * M_C + x_steel_Mo * M_Mo + x_steel_Si * M_Si;

  const Real M_O = 15.999;
  const Real M_spinel = M_Mn + 2 * M_Cr + 4 * M_O ;

  // Densities [g/cm^3]
  const Real rho_steel(7.67);
  const Real rho_spinel(4.93);

  //Average volume occupied by a molecule of MnCr2O4 [nm^3]
  const Real V_spinel = 1e21 * M_spinel / rho_spinel / Na;

  // Concentration of metal atoms in the steel [at/cm^3]
//  const Real C_steel = rho_steel * Na / M_steel;

  // C4 model parameters
  const Real nu(3600 * 1e13);           //Frequency of jump [/hr]
  const Real a(0.5);           // Length of jump [nm]
  const Real E_m_VMn(1.82);     // Mn vacancies energy of migration [eV]
  const Real E_m_h(1.7);       // Holes energy of migration [eV]

  // Interfaces concentrations [at/cm^3]
  const Real C_VMn_ox_g = 1e-4 * 1 / V_spinel ;
  const Real C_VCr_ox_g = 1e-4 * 2 / V_spinel ;
  const Real C_h_ox_g = 2 * C_VMn_ox_g + 3 * C_VCr_ox_g ;

  const Real C_VMn_ox_m = 1e-8;
  const Real C_h_ox_m = 1e-8;

  // Temperature fixed at 700C for now
  const Real temperature(973.15);

  // Mobilities
  const Real mu_VMn = 4 * pow(a, 2) * nu * (-2) / (Kb * temperature) *
                       exp(-E_m_VMn / (Kb * temperature));
  const Real mu_h = 4 * pow(a, 2) * nu * 1 / (Kb * temperature) *
                       exp(-E_m_h / (Kb * temperature));

  // Solving for eta
  const Real A = 8 * mu_VMn * C_VMn_ox_g - mu_h * C_h_ox_m;
  const Real B = mu_h * (C_h_ox_g - C_h_ox_m);
  const Real C = mu_h * C_h_ox_g - 8 * mu_VMn * C_VMn_ox_m;

  const Real eta = (- B - sqrt(pow(B, 2) - 4 * A * C)) / (2 * A);
  //std::cout<< "eta : " << eta << std::endl;

  const Real gamma_VMn = mu_VMn * Kb * temperature * log(eta) *
                          (C_VMn_ox_g - C_VMn_ox_m * pow(eta,-2)) / (1 - pow(eta,-2));

  const Real v_ox = - gamma_VMn * V_spinel / delta;

  std::cout<< "Velocity : " << v_ox << std::endl;

  return v_ox;
}
