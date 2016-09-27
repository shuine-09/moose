/*************************************************/
/*           DO NOT MODIFY THIS HEADER           */
/*                                               */
/*                     BISON                     */
/*                                               */
/*    (c) 2015 Battelle Energy Alliance, LLC     */
/*            ALL RIGHTS RESERVED                */
/*                                               */
/*   Prepared by Battelle Energy Alliance, LLC   */
/*     Under Contract No. DE-AC07-05ID14517      */
/*     With the U. S. Department of Energy       */
/*                                               */
/*  Written Vikash Mishra vkmphysics@hotmail.com */
/*      See COPYRIGHT for full restrictions      */
/*************************************************/

#ifndef MOXVELOCITY_H
#define MOXVELOCITY_H

#include "Material.h"

//Forward Declarations
class MOXVelocity;
class Function;

template<>
InputParameters validParams<MOXVelocity>();

class MOXVelocity : public Material
{

  public:
    MOXVelocity(const InputParameters & parameters);
    virtual ~MOXVelocity();

  protected:
    virtual void computeQpProperties();
    virtual void initQpStatefulProperties();

  private:
    //  bool _has_temp;
    const VariableValue & _temp;
    const VariableGradient & _grad_temp;
    const Real _limit;
    MaterialProperty<Real> & _pore_velocity;
};

#endif
