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

#ifndef VESICLEVOLUMECONSTRAINT_H
#define VESICLEVOLUMECONSTRAINT_H

#include "NonlocalKernel.h"
#include "Assembly.h"
#include "VesicleVolumeConstraintUserObject.h"

#include "libmesh/quadrature.h"

class VesicleVolumeConstraint : public NonlocalKernel
{
public:
  VesicleVolumeConstraint(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  /// new method for jacobian contributions corresponding to non-local dofs
  virtual Real computeQpNonlocalJacobian(dof_id_type dof_index);

  const VesicleVolumeConstraintUserObject & _shp;
  const Real & _shp_integral;
  const std::vector<Real> & _shp_jacobian;

  const std::vector<dof_id_type> & _var_dofs;

  Real _alpha_v;
};

template<>
InputParameters validParams<VesicleVolumeConstraint>();

#endif //VESICLEVOLUMECONSTRAINT_H
