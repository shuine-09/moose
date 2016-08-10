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
#ifndef XFEMWEIBULLMATERIAL_H
#define XFEMWEIBULLMATERIAL_H

#include "Material.h"

class XFEM;

//Forward Declarations
class XFEMWeibullMaterial;

template<>
InputParameters validParams<XFEMWeibullMaterial>();

/**
 * Empty material for use in simple applications that don't need material properties.
 */
class XFEMWeibullMaterial : public Material
{
public:
  XFEMWeibullMaterial(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties();
  virtual void computeQpProperties();

  MaterialProperty<Real> & _weibull_eta;
  MaterialProperty<Real> & _weibull_eta_old;

private:
  XFEM *_xfem;
  int _weibull_modulus;
  Real _specimen_volume;
  Real _specimen_material_property;
  Real _eta;
};

#endif //XFEMWEIBULLMATERIAL_H
