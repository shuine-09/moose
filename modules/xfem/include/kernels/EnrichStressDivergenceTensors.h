/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef ENRICHSTRESSDIVERGENCETENSORS_H
#define ENRICHSTRESSDIVERGENCETENSORS_H

#include "ALEKernel.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

//Forward Declarations
class EnrichStressDivergenceTensors;
class RankTwoTensor;
class RankFourTensor;

template<>
InputParameters validParams<EnrichStressDivergenceTensors>();

/**
 * EnrichStressDivergenceTensors mostly copies from StressDivergence.  There are small changes to use
 * RankFourTensor and RankTwoTensors instead of SymmElasticityTensors and SymmTensors.  This is done
 * to allow for more mathematical transparancy.
 */
class EnrichStressDivergenceTensors : public ALEKernel
{
public:
  EnrichStressDivergenceTensors(const InputParameters & parameters);

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
};

#endif //ENRICHSTRESSDIVERGENCETENSORS_H
