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

#ifndef PORECONTINUITY_H
#define PORECONTINUITY_H

#include "Kernel.h"

class PoreContinuity;

template<>
InputParameters validParams<PoreContinuity>();

/**
 * This kernel implements the spatial part of the pore transport equation,
 * dp/dt + div(ap) = 0
 */
class PoreContinuity : public Kernel
{
  public:
    PoreContinuity(const InputParameters & parameters);

    virtual ~PoreContinuity();

  protected:
    virtual Real computeQpResidual();
    virtual Real computeQpJacobian();
    virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  private:
    // Helper function called by convection_velocity().  Computes the
    // Temperature-squared denominator term.
    Real convection_scalar();

    // The convective velocity is given by convection_scalar() * grad(T)
    // To be called from computeQpResidual and computeQpJacobian terms.
    RealVectorValue convection_velocity();

    bool _has_temperature;
    const VariableValue & _T;
    const VariableGradient & _grad_T;
    bool _debug;
    //  const MaterialProperty<RealVectorValue> & _pore_velocity;
    const MaterialProperty<Real> & _pore_velocity;
    unsigned int _v_var;
    Real _alpha;
    Real _beta;

    // (constant) coefficient for diffusion of pores.
    Real _nu;
};


#endif /* PORECONTINUITY_H */
