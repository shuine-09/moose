//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "LevelSetTriMaterialBase.h"

// Forward Declarations

/**
 * Compute a Real material property for bi-materials problem (consisting of two different materials)
 * defined by a level set function
 *
 */
class LevelSetTriMaterialReal : public LevelSetTriMaterialBase
{
public:
  static InputParameters validParams();

  LevelSetTriMaterialReal(const InputParameters & parameters);

protected:
  virtual void assignQpPropertiesForLevelSetNegNeg();
  virtual void assignQpPropertiesForLevelSetPosNeg();
  virtual void assignQpPropertiesForLevelSetPosPos();

  /// Real Material properties for the three separate materials in the tri-material system
  std::vector<const MaterialProperty<Real> *> _bimaterial_material_prop;

  /// Global Real material property (switch tri-material diffusion coefficient based on level set values)
  MaterialProperty<Real> & _material_prop;
};
