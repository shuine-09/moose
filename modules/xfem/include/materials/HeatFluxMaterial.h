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
#ifndef HEATFLUXMATERIAL_H
#define HEATFLUXMATERIAL_H

#include "Material.h"
#include "MaterialProperty.h"

//Forward Declarations
class HeatFluxMaterial;

template<>
InputParameters validParams<HeatFluxMaterial>();

class HeatFluxMaterial : public Material
{
public:
  HeatFluxMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

  MaterialProperty<Real> & _mat_prop;
  Real _heat_transfer_coef;

  const VariableValue & _temp;
  const VariableValue & _temp_neighbor;
};

#endif //HEATFLUXMATERIAL_H
