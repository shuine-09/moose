/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef COMPUTELINEARELASTICSTRESSINTERFACE_H
#define COMPUTELINEARELASTICSTRESSINTERFACE_H

#include "Material.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "RotationTensor.h"
#include "DerivativeMaterialInterface.h"

class ComputeLinearElasticStressInterface : public DerivativeMaterialInterface<Material>
{
public:
  ComputeLinearElasticStressInterface(const InputParameters & parameters);
  virtual ~ComputeLinearElasticStressInterface() {}

protected:
  virtual void initQpStatefulProperties() override;
  virtual void resetQpProperties() override;
  virtual void computeQpProperties() override;

  /// Coupled displacement variables
  unsigned int _ndisp;
  std::vector<const VariableValue *> _disp;
  std::vector<const VariableGradient *> _grad_disp;

  std::string _base_name;

  MaterialProperty<RankTwoTensor> & _mechanical_strain;

  MaterialProperty<RankTwoTensor> & _total_strain;

  std::vector<MaterialPropertyName> _eigenstrain_names;
  std::vector<const MaterialProperty<RankTwoTensor> *> _eigenstrains;

  RankFourTensor _Cijkl;
  Real _poissons_ratio;
  Real _youngs_modulus;

  MaterialProperty<RankTwoTensor> & _stress;
};

#endif // COMPUTELINEARELASTICSTRESSINTERFACE_H
