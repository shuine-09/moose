//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "SideSetsGeneratorBase.h"

#include "libmesh/point.h"

// Forward declarations
class SideSetsFromAllElementFaces;

template <>
InputParameters validParams<SideSetsFromAllElementFaces>();

/**
 * Adds the faces on the boundary of given block
 * to the sidesets specified by "boundary"
 * Optionally, only adds faces that have a normal
 * equal to specified normal up to a tolerance
 */
class SideSetsFromAllElementFaces : public SideSetsGeneratorBase
{
public:
  SideSetsFromAllElementFaces(const InputParameters & parameters);

  std::unique_ptr<MeshBase> generate() override;

protected:
  std::unique_ptr<MeshBase> & _input;

  /// name of the sideset to which the interior faces will be added
  std::vector<BoundaryName> _boundary_names;
};
