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
#ifndef XFEMGAPFLUXMATERIAL_H
#define XFEMGAPFLUXMATERIAL_H

#include "Material.h"
#include "MaterialProperty.h"

//Forward Declarations
class XFEMGapFluxMaterial;

template<>
InputParameters validParams<XFEMGapFluxMaterial>();

class XFEMGapFluxMaterial : public Material
{
public:
  XFEMGapFluxMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

  MaterialProperty<Real> & _mat_prop;
  //MaterialProperty<Real> & _mat_prop_old;

  Real _heat_transfer_coef;
  Real _gap_tolerance;

  const VariableValue & _temp;
  const VariableValue & _temp_neighbor;

  const VariableValue & _disp_x;
  const VariableValue & _disp_x_neighbor;
  const VariableValue & _disp_y;
  const VariableValue & _disp_y_neighbor;
};

#endif //XFEMGAPFLUXMATERIAL_H
