//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef COMPUTEMULTIPHASELINEARELASTICPFFRACTURESTRESS_H
#define COMPUTEMULTIPHASELINEARELASTICPFFRACTURESTRESS_H

#include "ComputeIsotropicLinearElasticPFFractureStress.h"

class ComputeMultiPhaseLinearElasticPFFractureStress;

template <>
InputParameters validParams<ComputeMultiPhaseLinearElasticPFFractureStress>();

/**
 * Phase-field fracture
 * This class computes the stress and energy contribution to fracture
 * Small strain Anisotropic Elastic formulation
 * Stiffness matrix scaled for heterogeneous elasticity property
 */
class ComputeMultiPhaseLinearElasticPFFractureStress : public ComputeLinearElasticPFFractureStress
{
public:
  ComputeMultiPhaseLinearElasticPFFractureStress(const InputParameters & parameters);

protected:
  virtual void computeQpStress();

  /// switching function name list
  std::vector<MaterialPropertyName> _h_list;

  /// number of phases handled by this material
  unsigned int _n_phase;

  /// switching functions
  std::vector<const MaterialProperty<Real> *> _h_eta;

  // phase material properties
  std::vector<std::string> _phase_base;
  std::vector<const MaterialProperty<RankTwoTensor> *> _phase_stress;
  std::vector<const MaterialProperty<RankFourTensor> *> _dphase_stress_dstrain;

  // global material properties
  std::string _base_name;
};

#endif // COMPUTEMULTIPHASELINEARELASTICPFFRACTURESTRESS_H
