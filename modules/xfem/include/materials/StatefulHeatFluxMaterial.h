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
#ifndef STATEFULHEATFLUXMATERIAL_H
#define STATEFULHEATFLUXMATERIAL_H

#include "Material.h"
#include "MaterialProperty.h"

//Forward Declarations
class StatefulHeatFluxMaterial;

template<>
InputParameters validParams<StatefulHeatFluxMaterial>();

class StatefulHeatFluxMaterial : public Material
{
public:
  StatefulHeatFluxMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();
  virtual void initQpStatefulProperties();

  MaterialProperty<Real> & _mat_prop;
  MaterialProperty<Real> & _mat_prop_old;
  Real _heat_transfer_coef;

  const VariableValue & _temp;
  const VariableValue & _temp_neighbor;
};

#endif //STATEFULHEATFLUXMATERIAL_H
