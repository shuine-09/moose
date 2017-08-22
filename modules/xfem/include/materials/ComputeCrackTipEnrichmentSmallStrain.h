/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef COMPUTECRACKTIPENRICHMENTSMALLSTRAIN_H
#define COMPUTECRACKTIPENRICHMENTSMALLSTRAIN_H

#include "Material.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "RotationTensor.h"
#include "DerivativeMaterialInterface.h"
#include "Assembly.h"

/**
 * ComputeCrackTipEnrichmentSmallStrain calcualte the sum of standard strain and enrichement strain
 */
class ComputeCrackTipEnrichmentSmallStrain : public DerivativeMaterialInterface<Material>
{
public:
  ComputeCrackTipEnrichmentSmallStrain(const InputParameters & parameters);
  virtual ~ComputeCrackTipEnrichmentSmallStrain() {}

protected:
  virtual void initQpStatefulProperties();

  virtual void computeQpProperties();

  std::vector<Real> _enrich_disp;
  std::vector<RealVectorValue> _grad_enrich_disp;

  std::vector<std::vector<MooseVariable *>> _enrich_variable;

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

#endif // COMPUTECRACKTIPENRICHMENTSMALLSTRAIN_H
