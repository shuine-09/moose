/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CRACKTIPENRICHMENTSTRESSDIVERGENCETENSORS_H
#define CRACKTIPENRICHMENTSTRESSDIVERGENCETENSORS_H

#include "ALEKernel.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "CrackFrontDefinition.h"

// Forward Declarations
class CrackTipEnrichmentStressDivergenceTensors;
class RankTwoTensor;
class RankFourTensor;
class XFEM;

template <>
InputParameters validParams<CrackTipEnrichmentStressDivergenceTensors>();

/**
 * CrackTipEnrichmentStressDivergenceTensors implements the residual and jacobian for enrichement
 * displacement variables.
 *
 */
class CrackTipEnrichmentStressDivergenceTensors : public ALEKernel
{
public:
  CrackTipEnrichmentStressDivergenceTensors(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  std::string _base_name;

  const MaterialProperty<RankTwoTensor> & _stress;
  const MaterialProperty<RankFourTensor> & _Jacobian_mult;

  const MaterialProperty<RankTwoTensor> * _deformation_gradient;
  const MaterialProperty<RankTwoTensor> * _deformation_gradient_old;
  const MaterialProperty<RankTwoTensor> * _rotation_increment;

  const unsigned int _component;
  const unsigned int _enrichment_component;

  /// Coupled displacement variables
  unsigned int _ndisp;
  const std::vector<NonlinearVariableName> & _nl_vnames;
  std::vector<unsigned int> _enrich_disp_var;

private:
  MooseSharedPointer<XFEM> _xfem;
  const CrackFrontDefinition * _crack_front_definition;
};

#endif // CRACKTIPENRICHMENTSTRESSDIVERGENCETENSORS_H
