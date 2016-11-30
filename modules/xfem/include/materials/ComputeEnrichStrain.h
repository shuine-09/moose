/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef COMPUTEENRICHSTRAIN_H
#define COMPUTEENRICHSTRAIN_H

#include "Material.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "RotationTensor.h"
#include "DerivativeMaterialInterface.h"
#include "Assembly.h"

/**
 * ComputeEnrichStrain is the base class for strain tensors
 */
class ComputeEnrichStrain : public DerivativeMaterialInterface<Material>
{
public:
  ComputeEnrichStrain(const InputParameters & parameters);
  virtual ~ComputeEnrichStrain() {}

protected:
  virtual void initQpStatefulProperties();

  virtual void computeQpProperties();
  
  std::vector<Real> _enrich_disp;
  std::vector<RealVectorValue> _grad_enrich_disp;

  std::vector<std::vector<MooseVariable *> > _enrich_variable;
  
  MaterialProperty<RankTwoTensor> & _enrich_strain;

  /// the current shape functions
  const VariablePhiValue & _phi;

  /// gradient of the shape function
  const VariablePhiGradient & _grad_phi;

  /// Coupled displacement variables
  unsigned int _ndisp;
  std::vector<const VariableValue *> _disp;
  std::vector<const VariableGradient *> _grad_disp;

  std::string _base_name;

  MaterialProperty<RankTwoTensor> & _mechanical_strain;

  MaterialProperty<RankTwoTensor> & _total_strain;

  const MaterialProperty<RankTwoTensor> & _eigenstrain;

};

#endif //COMPUTEENRICHSTRAIN_H
