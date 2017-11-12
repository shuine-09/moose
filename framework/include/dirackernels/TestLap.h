/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#pragma once

// Moose Includes
#include "DiracKernel.h"

//Forward Declarations
class TestLap;

template<>
InputParameters validParams<TestLap>();

/**
 * TOOD
 */
class TestLap : public DiracKernel
{
public:
  static InputParameters validParams();
  TestLap(const InputParameters & parameters);

  virtual void addPoints() override;
  virtual Real computeQpResidual() override;

protected:
  const VariablePhiSecond & _second_phi;
  const VariableTestSecond & _second_test;
  const VariableSecond & _second_u;
};
