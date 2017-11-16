/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef LEVELSETMULTISTRESSMATERIAL_H
#define LEVELSETMULTISTRESSMATERIAL_H

#include "Material.h"

// Forward Declarations
class LevelSetMultiStressMaterial;
class RankTwoTensor;
class RankFourTensor;

template <>
InputParameters validParams<LevelSetMultiStressMaterial>();

/**
 * Construct a global strain from the phase strains in a manner that is consistent
 * with the construction of the global elastic energy by DerivativeMultiPhaseMaterial.
 */
class LevelSetMultiStressMaterial : public Material
{
public:
  LevelSetMultiStressMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

  const VariableValue & _ls;

  // phase material properties
  std::vector<std::string> _phase_base;
  std::vector<const MaterialProperty<RankTwoTensor> *> _phase_stress;
  std::vector<const MaterialProperty<RankFourTensor> *> _dphase_stress_dstrain;

  // global material properties
  std::string _base_name;
  MaterialProperty<RankTwoTensor> & _stress;
  MaterialProperty<RankFourTensor> & _dstress_dstrain;
};

#endif // LEVELSETMULTISTRESSMATERIAL_H
