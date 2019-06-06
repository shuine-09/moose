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
#ifndef WEIBULLMATERIAL_H
#define WEIBULLMATERIAL_H

#include "Material.h"

// Forward Declarations
class WeibullMaterial;

template <>
InputParameters validParams<WeibullMaterial>();

/**
 * Empty material for use in simple applications that don't need material properties.
 */
class WeibullMaterial : public Material
{
public:
  WeibullMaterial(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties();
  virtual void computeQpProperties();

  MaterialProperty<Real> & _weibull;
  const MaterialProperty<Real> & _weibull_old;

private:
  int _weibull_modulus;
  Real _specimen_volume;
  Real _specimen_material_property;
  Real _eta;
};

#endif // WEIBULLMATERIAL_H
