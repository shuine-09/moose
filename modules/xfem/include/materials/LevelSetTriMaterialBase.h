//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"

// Forward Declarations
class XFEM;

/**
 * Base class for switching between materials in a tri-material system where the interfaces are defined
 * by level set functions.
 */
class LevelSetTriMaterialBase : public Material
{
public:
  static InputParameters validParams();

  LevelSetTriMaterialBase(const InputParameters & parameters);

protected:
  virtual void computeProperties();
  virtual void computeQpProperties();

  /**
   * assign the material properties for the negative-negative level set region.
   */
  virtual void assignQpPropertiesForLevelSetNegNeg() = 0;

  /**
   * assign the material properties for the positive-negative level set region.
   */
  virtual void assignQpPropertiesForLevelSetPosNeg() = 0;

  /**
   * assign the material properties for the positive-negative level set region.
   */
  virtual void assignQpPropertiesForLevelSetPosPos() = 0;

  /// global material properties
  const std::string _base_name;

  /// Property name
  std::string _prop_name;

  /// shared pointer to XFEM
  std::shared_ptr<XFEM> _xfem;

  /// The variable number of the level set variable we are operating on
  const unsigned int _ls_var_1_number;

  /// The variable number of the level set variable we are operating on
  const unsigned int _ls_var_2_number;

  /// system reference
  const System & _system1;

  /// the subproblem solution vector
  const NumericVector<Number> * _solution1;

  /// system reference
  const System & _system2;

  /// the subproblem solution vector
  const NumericVector<Number> * _solution2;

  /// use the negative-negative level set region's material properties
  bool _use_neg_neg_property;

  /// use the positive-negative level set region's material properties
  bool _use_pos_neg_property;
};
