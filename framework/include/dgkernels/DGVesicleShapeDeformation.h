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

#ifndef DGVESICLESHAPEDEFORMATION_H
#define DGVESICLESHAPEDEFORMATION_H

#include "DGKernel.h"

// Forward Declarations
class DGVesicleShapeDeformation;

template <>
InputParameters validParams<DGVesicleShapeDeformation>();

/**
 * DG kernel for vesilce shape deformation
 */
class DGVesicleShapeDeformation : public DGKernel
{
public:
  DGVesicleShapeDeformation(const InputParameters & parameters);

protected:
  Real _eta;
  Real _epsilon;
  Real _h_elem;
  bool _rz;

  const VariablePhiSecond & _second_phi;
  const VariableTestSecond & _second_test;
  const VariableSecond & _second_u;
  const VariablePhiSecond & _second_phi_neighbor;
  const VariableTestSecond & _second_test_neighbor;
  const VariableSecond & _second_u_neighbor;

  virtual Real computeQpResidual(Moose::DGResidualType type);
  virtual Real computeQpJacobian(Moose::DGJacobianType type);
};

#endif
