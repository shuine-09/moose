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

#ifndef VESICLEAREACONSTRAINT_H
#define VESICLEAREACONSTRAINT_H

#include "NonlocalKernel.h"
#include "Assembly.h"
#include "VesicleAreaConstraintUserObject.h"

#include "libmesh/quadrature.h"

class VesicleAreaConstraint : public NonlocalKernel
{
public:
  VesicleAreaConstraint(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  /// new method for jacobian contributions corresponding to non-local dofs
  virtual Real computeQpNonlocalJacobian(dof_id_type dof_index);

  const VesicleAreaConstraintUserObject & _shp;
  const Real & _shp_integral;
  const std::vector<Real> & _shp_jacobian;

  const std::vector<dof_id_type> & _var_dofs;

  Real _alpha_a;
  Real _epsilon;

  const VariableSecond & _second_u;

};

template<>
InputParameters validParams<VesicleAreaConstraint>();

#endif //VESICLEAREACONSTRAINT_H
