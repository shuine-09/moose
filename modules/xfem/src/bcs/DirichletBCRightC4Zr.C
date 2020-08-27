//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DirichletBCRightC4Zr.h"

registerMooseObject("MooseApp", DirichletBCRightC4Zr);

defineLegacyParams(DirichletBCRightC4Zr);

/**
 * Boundary condition of a Dirichlet type
 *
 * Sets the value in the node
 *
 * Specific to the weak discontinuity equivalent of the C4 model for the
 * high-temperature corrosion of Zircaloy-4 (1000C to 1500C)
 * The variable used is the weak (atomic) oxygen concentration scaled by
 * the atomic concentration of Zr in the metal.
 */

InputParameters
DirichletBCRightC4Zr::validParams()
{
  InputParameters params = DirichletBCBase::validParams();
  //params.addRequiredParam<Real>("value", "Value of the BC");
  params.addParam<bool>("two_interfaces",true,"Boolean specifying if there are two interfaces"
                                "or not (if not, there's only one)");
  params.addParam<Real>("temperature",1473.15, "Temperature [K] in the cladding.");
  //params.addParam<Real>("temperature_ab",1473.15,"Temperature [K] at the alpha/beta interface");
  //params.declareControllable("value");
  params.addClassDescription("Imposes the right boundary condition for the C4 model of Zr."
                             "$u=g$, where u is the reduced weak discontinuity concentration"
                             " $g$ is a constant, depending of the number and value of interfacial jumps.");
  return params;
}

DirichletBCRightC4Zr::DirichletBCRightC4Zr(const InputParameters & parameters)
  : DirichletBCBase(parameters),
    _two_interfaces(getParam<bool>("two_interfaces")),
    _temperature(getParam<Real>("temperature"))
    //_temperature_ab(getParam<Real>("temperature_ab"))
{}

Real
DirichletBCRightC4Zr::computeQpValue()
{ const Real Zr_PBR = 1.55;
  const Real x_ox = 0.66666666667;
  //std::cout<< x_ox << std::endl;
  const Real C_ox = 3/Zr_PBR * x_ox;

  const Real x_ox_a = 0.667118123 - 1.10606e-5 * _temperature;
  const Real C_ox_a = 3/Zr_PBR * x_ox_a;

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
  const Real C_a_ox = x_a_ox/(1-x_a_ox);

  if (_two_interfaces)
  {
    Real x_b_a = (9.59e-3 * (_temperature - 1136) + 4.72e-6 * pow(_temperature - 1136,2) - 4.35e-9 * pow(_temperature - 1136,3)) * 1e-2;
    Real x_a_b = (45.86e-3 * (_temperature - 1136) - 44.77e-6 * pow(_temperature - 1136,2) + 17.40e-9 * pow(_temperature - 1136,3)) * 1e-2; //the original one, not the weak equivalent
    const Real C_a_b = x_a_b / (1 - x_a_b);
    const Real C_b_a = x_b_a / (1 - x_b_a);

    _value = C_ox - (C_ox_a - C_a_ox) - (C_a_b - C_b_a);

    //std::cout << "C_ox : " << C_ox << std::endl;
    //std::cout << "C_ox_a : " << C_ox_a << std::endl;
    //std::cout << "C_a_ox : " << C_a_ox << std::endl;
    //std::cout << "C_a_b : " << C_a_b << std::endl;
    //std::cout << "C_b_a : " << C_b_a << std::endl;
  }
  else
  {
    _value = C_ox - (C_ox_a - C_a_ox);
  }
  return _value;
}
