/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef MAXIMUMNORMALSEPARATION3D_H
#define MAXIMUMNORMALSEPARATION3D_H

#include "InterfaceMaterial.h"

class MaximumNormalSeparation3D;

template <>
InputParameters validParams<MaximumNormalSeparation3D>();

/**
 *
 */
class MaximumNormalSeparation3D : public InterfaceMaterial
{
public:
  MaximumNormalSeparation3D(const InputParameters & parameters);

protected:
  virtual void resetQpProperties() override;
  virtual void computeQpProperties() override;
  virtual void initQpStatefulProperties() override;

  const std::string _base_name;
  MaterialProperty<Real> & _max_normal_separation;
  const MaterialProperty<Real> & _max_normal_separation_old;

  const VariableValue & _disp_x;
  const VariableValue & _disp_x_neighbor;
  const VariableValue & _disp_y;
  const VariableValue & _disp_y_neighbor;
  const VariableValue & _disp_z;
  const VariableValue & _disp_z_neighbor;
};

#endif // MAXIMUMNORMALSEPARATION3D_H
