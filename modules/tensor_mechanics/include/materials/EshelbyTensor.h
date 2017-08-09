/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef ESHELBYTENSOR_H
#define ESHELBYTENSOR_H

#include "Material.h"

class RankTwoTensor;

/**
 * EshelbyTensor defines a strain increment and rotation increment, for finite strains.
 */
class EshelbyTensor : public Material
{
public:
  EshelbyTensor(const InputParameters & parameters);

  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

protected:
  std::string _base_name;

  MaterialProperty<Real> & _sed;
  const MaterialProperty<Real> & _sed_old;
  MaterialProperty<RankTwoTensor> & _eshelby_tensor;
  const MaterialProperty<RankTwoTensor> & _stress;
  const MaterialProperty<RankTwoTensor> & _stress_old;
  const MaterialProperty<RankTwoTensor> & _strain_increment;

  std::vector<const VariableGradient *> _grad_disp;

private:
};

#endif // ESHELBYTENSOR_H
