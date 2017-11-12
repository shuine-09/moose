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


#include "DGKernel.h"

//Forward Declarations
class DGLinearStommelMunk;

/**
 * DG kernel for linear stommel munk
 */
class DGLinearStommelMunk : public DGKernel
{
public:
  DGLinearStommelMunk(const InputParameters & parameters);
  static InputParameters validParams();

protected:
  Real _eps_m;
  Real _eps_s;
  Real _eta;

  const VariablePhiSecond & _second_phi;
  const VariableTestSecond & _second_test;
  const VariableSecond & _second_u;
  const VariablePhiSecond & _second_phi_neighbor;
  const VariableTestSecond & _second_test_neighbor;
  const VariableSecond & _second_u_neighbor;

  virtual Real computeQpResidual(Moose::DGResidualType type) override;
  virtual Real computeQpJacobian(Moose::DGJacobianType type) override;

};
