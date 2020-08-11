//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "InitialCondition.h"

#include <string>

// Forward Declarations
class C4ZrICConst;
class InputParameters;

template <typename T>
InputParameters validParams();

template <>
InputParameters validParams<C4ZrICConst>();

/**
 * Defines a boundary condition that forces the value to be a user specified
 * function at the boundary.
 */
class C4ZrICConst : public InitialCondition
{
public:
  static InputParameters validParams();

  C4ZrICConst(const InputParameters & parameters);

protected:
  /**
   * Evaluate the function at the current quadrature point and time step.
   */
  //Real f();

  /**
   * The value of the variable at a point.
   */
  virtual Real value(const Point & p) override;

  /**
   * The value of the gradient at a point.
   */
  //virtual RealGradient gradient(const Point & p) override;

  Real _temperature;

//  Real _x_b_break = 531.0 ;

  Real _x_a_b = 572.2 ;

//  Real _x_a_break = 580.8;

  Real _x_ox_a = 590.0 ;

  Real _C_b = 0.0074/(1-0.0074) ;

  Real _C_b_a = 0.0373 ;

  Real _C_a_b = 0.1241 ;

  Real _C_a_ox = 0.4241 ;

  Real _C_ox_a = 1.2597;

//  Real _C_ox_a_weak = 0.3373;

  Real _C_ox = 1.2903;

  Real _C_ox_weak = 0.3679;

//  Real _C_a_break = 0.0907;

//  Real _grad_ba = 723.17 * 1e-6;

//  Real _grad_ab = 6249.18 * 1e-6;

//  Real _grad_aox = 26620.03 * 1e-6;

//  Real _grad_ox = (_C_ox-_C_ox_a)/(600-_x_ox_a);
};
